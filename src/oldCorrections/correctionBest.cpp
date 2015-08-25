#include <random>
#include <array>
#include <functional>

#include "Timer.h"
//#include "trie.h"
#include "classIO.h"

//#define DEBUG 1
#define MAX_DEPTH 5
#define MAX_KMERS(x)  (100 - x + 1)

//NEW IN THIS VERSION
//4. Removing loop unroling from correctReads() and do all the messy stuf in dfs()
//3. Multi-threading
//2. Estimating depth based on the distance between c-sets
//1. searchUpdate() only called at the end, after correcting the k-mers 

unsigned char c[8]={ 0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80};
unsigned char r[8]={ 0xFE, 0xFD, 0xFB, 0xF7, 0xEF, 0xDF, 0xBF, 0x7F};

unsigned int hashSizes[SIZES]={30000001, 40000003, 50000017, 60000011,70000027, 80000023, 90000049, 100000567, 150000001, 200000477, 250000013, 300000641, 350000041, 400000651, 450000007, 500000461, 550000001, 600000607, 650000011, 700000747, 750000007, 800000533, 850000013, 900000637, 950000017, 1000000531, 1250000027, 1500000713, 1750000027, 2000000579, 2147483647}; 

using namespace std;
using namespace general;
using namespace IOFn;

//namespace correction {

unsigned short cores;
unsigned int kmerSize;

//all pools data
IOF<mybitsetx>* ioh;
Trie* trie;
HashTable<mysig>* pools2bacs;
	
//per pool data
//string readFname;
ofstream out;
	
unsigned int numReads;
struct readFasta* reads;
	
bool correct(char* read, unsigned short direction, int crtKmer, int pos, int kmerPos, mybitsetx c, vector<mybitsetx>& kmersInHash, const vector<unsigned short>& orientation, vector<cset>& csets, unsigned short crtCset, unsigned int& numSearch, unsigned int& numHit) {

bool correct=false;
//mybitsetx e=kmers[crtKmer];
//mybitsetx e=kmersInHash[crtKmer];
unsigned short bacs[MAX_BACS];
unsigned short pools[POOLS];
unsigned short numBacs, numPools; 

unsigned short orient=orientation[crtKmer];

mybitsetx g;
mybitsetx o(c);
c.invert();
if (c<o)
	g=c;
else
	g=o;

//DEBUG
//cout << "[correct]Searching " << g; 
//#pragma omp critical
//{
if (!ioh->searchCopy(g)) {
//}
#if defined DEBUG
	cout << "[correct]kmer not found in the hashtable!" << endl;
#endif	
	return false;
};	

numPools=g.getPools(pools);
if (numPools<LOW) {
#if defined DEBUG
	cout << "[correct]No pools found " << numPools << endl;
#endif	
	return false;	
};

if (numPools>HIGH)
	correct=true;

numBacs=0;
if (numPools>=LOW&&numPools<=HIGH) {
	mysig sig(pools, numPools);
	if (!pools2bacs->searchCopy(sig)) { 
		trie->search(pools, numPools, bacs, numBacs);
		
		numSearch++;
		
		sig.setBacs(bacs, numBacs);
		#pragma omp critical (insert)
		{
		pools2bacs->insert(sig);
		}
	} else {
		numHit++;
	}
	numBacs=sig.getNumBacs();	
	for (unsigned short i=0;i<numBacs;i++)
		bacs[i]=sig.getBac(i);
	
	if (numBacs==0)	{
#if defined DEBUG
		cout << "[correct]No bacs found " << numBacs << endl;
		cout << "[correct]numPools " << numPools << ": ";
		for (unsigned short i=0;i<numPools;i++)
			cout << pools[i] << " ";
		cout << endl;	
#endif		
		return false;
	}
	
	bool hasBacs=false, found=false;
	//check that bacs found agree with other cset bacs
	for (auto it=csets.begin();it!=csets.end();it++) {
		if (it->bacs.size()==0)
			continue;
			
		hasBacs=true;
		
		//found=false;
		for (unsigned short i=0;i<numBacs;i++) {//does this kmer belong to crt cset?
			if (it->bacs.find(bacs[i])!=it->bacs.end()) {
				found=true;
				break;
			}
		}
		
		if (found)
			break;
	}
	
	if (!hasBacs||found) {
		correct=true;
	} else {
#if defined DEBUG
		cout << "[correct]Bacs do not match: ";
		for (unsigned short i=0;i<numBacs;i++) 
			cout << bacs[i] << ",";
		cout << endl;
		return false;
#endif		
	}

	if (found) {
		set<unsigned short> newBacs;
		set<unsigned short> crtBacs=csets[crtCset].bacs;
		for (unsigned short i=0;i<numBacs;i++)
			if (crtBacs.find(bacs[i])!=crtBacs.end()) 
				newBacs.insert(bacs[i]);
		csets[crtCset].bacs=newBacs;		
	} else {
		csets[crtCset].bacs.insert(bacs, bacs + numBacs);
	}
	
}//if (numPools>=LOW&&numPools<=HIGH) {

#if defined DEBUG
cout << "[correct]Found " << g; 
cout << "[correct]numPools " << numPools << ": ";
for (unsigned short i=0;i<numPools;i++) 
	cout << pools[i] << ",";
cout << endl;
cout << "[correct]numBacs " << numBacs << ": ";
for (unsigned short i=0;i<numBacs;i++) 
	cout << bacs[i] << ",";
cout << endl;
#endif

if (correct) {
	
	//change base in read
	char newBase;
#if defined DEBUG
	char oldBase=read[pos];
#endif
	if (orient)
		newBase=o.getBase(kmerPos);
	else 
		newBase=mybitsetx::rev(o.getBase(kmerPos)); 
	
	read[pos]=newBase;
#if defined DEBUG
	cout << "[correct]read pos " << pos << " old base " << oldBase << " new base " << newBase << endl;
#endif
	
	//change original k-mer so that we can reverse pool counts in hash table if this correction turns out to be wrong and we have to backtrack
	kmersInHash[crtKmer]=g;
	
	if (direction==0) 
		csets[crtCset].begin--;
	else	
		csets[crtCset].end++;
	
#if defined DEBUG
	cout << "[correct]crt cset begin " << csets[crtCset].begin << ", crt cset end " << csets[crtCset].end << endl << endl; 
#endif
	
	return true;

}//if (correct)
	
return false;
};

void dfs(char* read, char* original, const vector<mybitsetx>& kmers, vector<mybitsetx>& kmersInHash, const vector<unsigned short>& orientation, vector<cset> csets, unsigned int crt, unsigned short depth, vector<unsigned short> posChanged, unsigned short& minCorrections, unsigned short& numCorrections,  unsigned int& numSearch, unsigned int& numHit) {
	
unsigned int num=strlen(read);	

numCorrections=posChanged.size();

mybitsetx c, e;
unsigned short orient=0;
unsigned short direction=0;
int crtKmer=0, pos=0, kmerPos=0;

int prevECset=-1;
int nextBCset=num - kmerSize + 1;

unsigned int crtCset=crt;

while (true) {
	
	prevECset=-1;
	if (crtCset>0)
		prevECset=csets[crtCset - 1].end;
	nextBCset=num - kmerSize + 1;
	if (crtCset + 1 < csets.size())
		nextBCset=csets[crtCset + 1].begin;
	
	while (((csets[crtCset].begin - 1)==prevECset)&&((csets[crtCset].end + 1)==nextBCset)) {
		
		crtCset++;
		if (crtCset==csets.size())
			break;
		
		prevECset=csets[crtCset-1].end;
			
		nextBCset=num - kmerSize + 1;
		if (crtCset + 1 < csets.size())
			nextBCset=csets[crtCset + 1].begin;
				
	}//end while			
	
	if (crtCset==csets.size())
		break;
	
	//update direction
	if ((csets[crtCset].begin - 1)!=prevECset) 
		direction=0;
	else if ((csets[crtCset].end + 1)!=nextBCset) 
		direction=1;
		
	
	//determine pos in read&kmer which need to be changed next
	if (direction==0) {
		crtKmer=csets[crtCset].begin - 1;
		assert(crtKmer>=0);
		//cout << "direction 0 " << ", crtKmer " << crtKmer << endl; 
		pos=crtKmer;
		orient=orientation[crtKmer];
		if (orient)
			kmerPos=0;
		else
			kmerPos=kmerSize - 1;
	} else { 
		crtKmer=csets[crtCset].end + 1;
		assert(crtKmer<num - kmerSize + 1);
		//cout << "direction 1 " << ", crtKmer " << crtKmer << endl; 
		pos=crtKmer + kmerSize - 1;
		orient=orientation[crtKmer];
		if (orient)
			kmerPos=kmerSize - 1;
		else
			kmerPos=0;
	}	
	
	//make a copy of the original k-mer
	if (!(crtKmer>=0)||!(crtKmer<num - kmerSize + 1)) {
		cout << "crtKmer " << crtKmer << endl; exit(1);
	}
	c=kmers[crtKmer];
	
#if defined DEBUG
cout << "[dfs]crtKmer " << crtKmer << endl;
cout << "[dfs]num corrections so far " << numCorrections << endl;
if (posChanged.size()>0)
	cout << "[dfs]Updating kmer with existing read changes" << endl;
#endif

	int p, kp;
#if defined DEBUG		
	char oldBase;
#endif			
	//update crt kmer with possible base changes from previous fixed kmers
	for (auto it2=posChanged.begin();it2!=posChanged.end();it2++) {
		p=(*it2);
		//assert(p>=0&&p<(int)num);
		if ((crtKmer<=p)&&(p<=(crtKmer + int(kmerSize) - 1))) {
			if (orient) {
				kp=p - crtKmer;
#if defined DEBUG		
				oldBase=c.getBase(kp);
#endif		
				c.setBase(kp, read[p]);
			} else {
				kp=kmerSize - 1 - (p - crtKmer);
#if defined DEBUG		
				oldBase=c.getBase(kp);
#endif		
				c.setBase(kp, mybitset::rev(read[p]));
			}	
#if defined DEBUG
cout<<"\t[dfs](update)crtKmer "<<crtKmer<<" read pos prev changed "<<p<< " kmerPos updated "<<kp<<" old base "<<oldBase<<" new base "<<c.getBase(kp)<<endl<<endl;   
#endif
		}//end if
	}//end for

#if defined DEBUG
cout << "[dfs]Direction " << direction << " read pos to change " << pos << " orientation " << orient << " kmerPos to change " << kmerPos << endl;
cout << "[dfs]read ";
for (int j=crtKmer;j<=(crtKmer+int(kmerSize)-1);j++)
	cout << read[j];
cout << endl;	
cout << "[dfs]kmer " << c;
#endif	
		
#if defined DEBUG
cout << "[dfs]Trying same k-mer " << c;
#endif	
	bool corrected=correct(read, direction, crtKmer, pos, kmerPos, c, kmersInHash, orientation, csets, crtCset, numSearch, numHit);
		
	if (!corrected)
		break;
		
}//end while (true)		
	

if (crtCset==csets.size()) {
//	if (numCorrections==depth) {

#if defined DEBUG
cout<<"[dfs]Found SOLUTION PATH with cost " << numCorrections <<endl<<endl;
#endif
		minCorrections=numCorrections;
//	} 	
	return;
}//end if
	
if ((numCorrections + 1)>depth) {
#if defined DEBUG
cout<<"[dfs]Exceded MAX DEPTH. Returning."<<endl<<endl;
#endif
	return;
}

//recurse
vector<char> toTry;
//toTry.push_back(c.getBase(kmerPos));
c.flip(2*kmerPos);
toTry.push_back(c.getBase(kmerPos));
c.flip(2*kmerPos);
c.flip(2*kmerPos + 1);
toTry.push_back(c.getBase(kmerPos));
c.flip(2*kmerPos);
toTry.push_back(c.getBase(kmerPos));

for (unsigned int i=0;i<toTry.size();i++) { 
	
	c.setBase(kmerPos, toTry[i]);
#if defined DEBUG
cout << "[dfs]crtKmer " << crtKmer << endl;
cout << "[dfs]Trying " << c;
#endif	
	
	bool corrected=correct(read, direction, crtKmer, pos, kmerPos, c, kmersInHash, orientation, csets, crtCset,numSearch,numHit);
	
	if (corrected) {
		//record read pos which was succesfully changed
		//if (read[pos]!=original[pos])
		posChanged.push_back(pos);
		//dfs(read, original, kmers, kmersInHash, orientation, csets, crtCset, depth, posChanged, minCorrections, numCorrections, numSearch, numHit);
		if (minCorrections==depth)
			break;
	}		
	
}//end for
	
};//end method definition 


void correctReads(unsigned int pool, unsigned int& numSearch, unsigned int& numHit) {

//#if defined DEBUG
struct kmerInfo info[MAX_KMERS(kmerSize)];
//#endif

unsigned short orient;

vector<mybitsetx> kmers;
vector<mybitsetx> kmersInHash;
vector<unsigned short> orientation;

unsigned short numPools;
unsigned short pools[POOLS];
unsigned short numBacs;
unsigned short bacs[MAX_BACS];

omp_set_num_threads(cores);
printf("Number of cores is %d\n", cores);	

unsigned int noCset=0;
unsigned int totalCorrect=0;
unsigned int nonRepReads=0;
unsigned int canCorrect=0;
unsigned int cannotCorrect=0;
unsigned int totalChanged=0;
unsigned int totalUnchanged=0;

//map<unsigned int, ostringstream*> outs;
//for (unsigned int i=0;i<cores;i++) {
//	outs[i]=new ostringstream();
//}
#pragma omp parallel default(none) shared(noCset, totalCorrect, nonRepReads, canCorrect, cannotCorrect, totalChanged,totalUnchanged,cout,kmerSize,pool,ioh,trie,pools2bacs,reads,numReads,out,numSearch,numHit) private(numPools,pools,numBacs,bacs,orient,kmers,kmersInHash,orientation,info) 
{
ostringstream os;
//unsigned int id = omp_get_thread_num();
//#pragma omp critical
//outs[id]=s;
#pragma omp for schedule(static)
for (unsigned int r=0;r<numReads;r++) { 
	
	char* header = reads[r].header;
	//cout << "header " << header << endl;
	char* read = reads[r].read;
	unsigned int num = strlen(read);
	char* original=(char*)malloc((num+1)*sizeof(char));
	strcpy(original, read);

	//count the number of errors the read contains
	char* errs=strchr(header, ':');
	unsigned int numErrors=0;
	while(*errs) if (*errs++ == '_') numErrors++;
	
	kmers.clear();
	kmersInHash.clear();
	orientation.clear();
	//corrections.clear();
	
	unsigned int crtCset;
	vector<cset> csets;
	vector<cset> _csets;

	bool conflicts=false;
	
	class mybitsetx b,d,e,g;
	for (unsigned int k=0; k<num-kmerSize+1; k++) {
		int z=0;
		int j=kmerSize*2-1;
		//bool findN=false;
		unsigned int c;
		for (c=0; c<kmerSize; c++) {
			switch (read[k+c]) {
			case 'A':
			b.set(z,0);
			z++;
			b.set(z,0);
			z++;
			d.set(j,1);
			j--;
			d.set(j,1);
			j--;
			break;
			case 'C':
			b.set(z,0);
			z++;
			b.set(z,1);
			z++;
			d.set(j,0);
			j--;
			d.set(j,1);
			j--;
			break;
			case 'G':
			b.set(z,1);
			z++;
			b.set(z,0);
			z++;
			d.set(j,1);
			j--;
			d.set(j,0);
			j--;
			break;
			case 'T':
			b.set(z,1);
			z++;
			b.set(z,1);
			z++;
			d.set(j,0);
			j--;
			d.set(j,0);
			j--;
			break;
			case 'N':
			b.set(z,0);
			z++;
			b.set(z,0);
			z++;
			d.set(j,1);
			j--;
			d.set(j,1);
			j--;
			read[k+c]='A';
			break;
			//default:
			//findN=true;
			//c=kmerSize;
			}
	}//end for c		
	
	//if (findN) 
	//	continue;
	
	if (d<b) 
		{e=d;orient=0;}
	else
		{e=b;orient=1;}
	
	//#pragma omp critical
	//{
	ioh->searchCopy(e);
	//}
	kmers.push_back(e);
	orientation.push_back(orient);
	
	//get positive pools and their number for this k-mer
	numPools=e.getPools(pools);

	numBacs=0;
	if (numPools>=LOW && numPools<=HIGH) { //get bacs for this kmer 
		mysig sig(pools, numPools);
		if (!pools2bacs->searchCopy(sig)) { 
			trie->search(pools, numPools, bacs, numBacs);
			
			numSearch++;
			
			sig.setBacs(bacs, numBacs);
			#pragma omp critical (insert)
			{
			pools2bacs->insert(sig);
			}
		} else {
			numHit++;
		}
		numBacs=sig.getNumBacs();	
		for (unsigned int i=0;i<numBacs;i++)
			bacs[i]=sig.getBac(i);
	}

#if defined DEBUG
	info[k].pos=k;
	info[k].numPools=numPools;
	info[k].numBacs=numBacs;
	for (unsigned int i=0;i<numBacs;i++)
		info[k].bacs[i]=bacs[i];
	for (unsigned int i=0;i<numPools;i++)
		info[k].pools[i]=pools[i];
#endif

	if (numPools<LOW||(numPools>=LOW&&numPools<=HIGH&&numBacs==0))
		continue;
	
	crtCset=_csets.size() - 1;

	if (numPools>HIGH) { 
		if (_csets.size()>0&&_csets[crtCset].bacs.size()==0&&k==(unsigned int)(_csets[crtCset].end + 1)) { 
			_csets[crtCset].end=k;
		} else {//create new cset	
			cset crt;
			crt.begin=k;
			crt.end=k;
			_csets.push_back(crt);	
#if defined DEBUG
			cout << "\t[correctReads]new cset @ k " << k << ", numPools " << numPools << endl;
#endif		
		}	
		
		continue;
	}	
	
	bool hasBacs=false, found=false;
	for (auto it=_csets.begin();it!=_csets.end();it++) {
		
		if (it->bacs.size()==0)
			continue;
		
		hasBacs=true;
		
		for (unsigned i=0;i<numBacs;i++) {//does this kmer belong to crt cset?
			if (it->bacs.find(bacs[i])!=it->bacs.end()) {
				found=true;
				break;
			}
		}

		if (found)  
			break;
	}

	//if (!hasBacs)
	//	found=true;
	
	if (found&&(_csets[crtCset].bacs.size()>0)&&k==(unsigned int)(_csets[crtCset].end + 1)) {//same cset	
		_csets[crtCset].end=k;
		//cset only contains the common bacs among all k-mers
		set<unsigned short> newBacs;
		set<unsigned short> crtBacs=_csets[crtCset].bacs;
		for (unsigned int i=0;i<numBacs;i++)
			if (crtBacs.find(bacs[i])!=crtBacs.end()) 
				newBacs.insert(bacs[i]);
		_csets[crtCset].bacs=newBacs;		
#if defined DEBUG
	//cout << "\t[correctReads]crt cset begin " << _csets[crtCset].begin << ", crt cset end " << _csets[crtCset].end << endl; 
#endif
	} else {//create new cset	
		cset crt;
		crt.begin=k;
		crt.end=k;
		crt.bacs.insert(bacs, bacs + numBacs);	
		_csets.push_back(crt);	
#if defined DEBUG
		cout << "\t[correctReads]new cset @ k " << k << endl;
#endif		
	}	
	
	if (!found&&hasBacs) 
		conflicts=true;

	}//end for k

#if defined DEBUG
cout << endl << header << ", numErrors "<< numErrors << ", numKmers " << (num - kmerSize + 1) << endl << read << endl;
for (unsigned int k=0;k<num-kmerSize+1;k++) {
	cout << "kmer " << info[k].pos << " numPools " << info[k].numPools << " numBacs " << info[k].numBacs << ", bacs ";
	for (unsigned int i=0;i<info[k].numBacs;i++) 
		cout << info[k].bacs[i] << ",";
	//if (info[k].numPools>=3 && info[k].numPools<=11) {
		//cout << " pools ";
		//for (unsigned int i=0;i<info[k].numPools;i++) 
		//	cout << info[k].pools[i] << ",";
	//}
	cout << endl;	
}
cout << header << " numErrors "<< numErrors << endl << endl;
#endif

vector<cset>::iterator prev=_csets.begin();
vector<cset>::iterator crt=_csets.begin();
vector<cset>::iterator next=_csets.begin();

bool nonRep=false;
unsigned short depth=1;
while (true) {
	//FIXME
	//if (((next!=_csets.end())&&(crt->end + 1 + (int)kmerSize < next->begin)) || ((next==_csets.end())&&(crt->end + (int)kmerSize < (int)num)))
	//if ((_csets[0].begin >0) || ((next!=_csets.end())&&(crt->end + 1 < next->begin)) || ((next==_csets.end())&&(crt->end + 1 < (int)(num - kmerSize + 1))))
	//	depth++;
	
	prev=crt;
	crt=next;
	if (crt==_csets.end())
		break;
	next++;
	
	int prevEnd=-1;
	if (prev!=_csets.begin())
		prevEnd=prev->end;
	int nextBegin=num-kmerSize+1;
	if (next!=_csets.end())
		nextBegin=next->begin;
	if ((crt->end - crt->begin + 1 <3)&&(crt->begin > prevEnd + 1)&&(crt->end + 1 < nextBegin))
		continue;

	//if (conflicts) {
	//if (conflicts&&crt->bacs.size()>0) {
	//	crt->end=crt->begin - 1;
	//	crt->bacs.clear();
	//}	
	
	//DEBUG
	if (crt->bacs.size()>0)
		nonRep=true;
	
	csets.push_back(*crt);

}//end while	

if (csets.size()==0) {
	#pragma omp atomic
	noCset++;
}	

if (nonRep) {
#pragma omp atomic
nonRepReads++;
}

kmersInHash=kmers;

//use csets to fix errors by moving left and right

#if defined DEBUG
cout << "CORRECTION" << endl;
#endif
	crtCset=0;
	unsigned int oldCset=crtCset;
	unsigned int nextCset=crtCset + 1;

	int prevECset=-1, nextBCset=num-kmerSize+1;
	while (crtCset!=csets.size()) {
		
		//update crt cset ends & neighbors	
		if (crtCset!=oldCset) 
			prevECset=csets[oldCset].end;
		
		nextBCset=num - kmerSize + 1;
		if (nextCset!=csets.size()) 
			nextBCset=csets[nextCset].begin;
				
		if (((csets[crtCset].begin - 1)==prevECset)&&((csets[crtCset].end + 1)==nextBCset)) {
			oldCset=crtCset;
			crtCset=nextCset;
			//if (nextCset!=csets.size())
				nextCset++;
			continue;	
		}

		break;

	}//end while
				
	//for correct reads don't do anything
	if (crtCset==csets.size()) { 
#if defined DEBUG
		cout << "[correctReads]Read already correct" << endl;
		cout <<"[correctReads]" << header << ":" << 1 << endl << endl;
#endif
		//os << header << ":" << 1 << endl << read << endl;
		//os << header << endl << read << endl;
		
		if (csets.size()>0) {
			#pragma omp atomic
			totalCorrect++;
		}	

		continue;
	};	
	
	//string correction;
	unsigned short numCorrections=0;
	unsigned short minCorrections=100;
	vector<unsigned short> posChanged;
	
	for (unsigned int crtCset=0;crtCset<csets.size();crtCset++) {
		depth=1;
		//FIXME
		vector<cset> crtCsets;
		if (!conflicts) 
			crtCsets=csets;
		else
			crtCsets.push_back(csets[crtCset]);
		while (depth<=MAX_DEPTH) {
#if defined DEBUG
cout << endl << "[correctReads]Crt depth is " << depth << endl << endl;
#endif	

			dfs(read, original, kmers, kmersInHash, orientation, crtCsets, crtCset, depth, posChanged, minCorrections, numCorrections, numSearch, numHit);
			
			if (minCorrections==depth)
				break;
	
			depth++;	
		}//end while 
		
		if (minCorrections==depth)
			break;
		
		if (!conflicts) 
			break;
		
	}//end for
	
	if (minCorrections==depth) {

#if defined DEBUG
cout << "[correctReads]Cost of best path is  " << minCorrections << endl;
cout << "[correctReads]" << header << ":" << 1 << endl << read << endl << endl;
#endif
		//os << header << ":" << 1 << endl << read << endl;
		//os << header << endl << read << endl;
		
		#pragma omp atomic
		canCorrect++;
	} else {

#if defined DEBUG
cout << "[correctReads]Unable to correct " << endl;
cout << "[correctReads]" << header << ":" << 0 << endl << endl;
#endif
		//os << header << ":" << 0 << endl << read << endl;
		//os << header << endl << read << endl;
		
		#pragma omp atomic
		cannotCorrect++;
		if (numCorrections>0) {
			#pragma omp atomic 
			totalChanged++;
		} else {
			os << header << endl << read << endl;
			#pragma omp atomic 
			totalUnchanged++;
		}
		continue;
	}//end if
	

//update k-mer counts in the hashtable
for (unsigned int k=0; k<num-kmerSize+1; k++) {
	e=kmers[k];
	g=kmersInHash[k];
	if (g!=e) {
		//#pragma omp critical (update)
		//{
		ioh->searchUpdatePool(e, pool, -1);
		ioh->searchUpdatePool(g, pool, 1);
		//}
	}	
}

}//end for r

#pragma omp critical 
out << os.str();

}//#pragma omp parallel 

cout << "[correctReads]noCset " << noCset << endl;
cout << "[correctReads]Total correct " << totalCorrect << endl;
cout << "[correctReads]Can correct " << canCorrect << endl;
cout << "[correctReads]Cannot correct " << cannotCorrect << endl;
cout << "[correctReads]Number of reads changed " << totalChanged << endl;
cout << "[correctReads]Number of reads unchanged " << totalUnchanged << endl;
cout << "[correctReads]Reads non-repetitive " << nonRepReads << endl;

};

int main(int argc, char* argv[]) {

if (argc!=7)
	{
	std::cerr<<"\n\nUSE: corr <path/ReadFile> <path/MapPoolsToBacsFile> <path/OutputFile> <k-mer_size> <input_file_number> <thread>\n\n";
	exit(EXIT_FAILURE);
	}

string readFname(argv[1]);
kmerSize=atoi(argv[4]);
cores=atoi(argv[6]); 
unsigned int files=atoi(argv[5]); 

char delimC[] = ".";
Parser parser;
parser.update(delimC, readFname);

cout << "MAX_KMERS is " << MAX_KMERS(kmerSize) << endl;

Timer* pTimer;

ioh=new IOF<mybitsetx>(argv[1], argv[2], argv[3], kmerSize, files);
ioh->ReadB();
//ioh.CheckCorrect();

//load BAC sigs into trie
unsigned short* bacPools=(unsigned short*)malloc(sizeof(unsigned short)*BACS*LAYERS);
readBacMappings(bacPools);

trie=new Trie();
for (unsigned short i=0;i<BACS;i++)
	trie->insert(&bacPools[i*LAYERS], 7, i+1);

pools2bacs=new HashTable<mysig>();
pools2bacs->allocHash(1);

reads=(struct readFasta*)malloc(MAX_POOL_SIZE*sizeof(struct readFasta));
if (!reads) {
	perror("Error allocating the reads!");
	exit(1);
}
for (unsigned int r=0; r<MAX_POOL_SIZE; r++) {
	reads[r].header=NULL;
	reads[r].read=NULL;
}

for (unsigned int i=60;i<61;i++) {
	
	numReads = ioh->readReads(i, reads);
	printf("Number of reads read is %d\n", numReads);  
	
	unsigned int numHit=0;
	unsigned int numSearch=0;
	
	//output files for corrected and uncorrected reads
	ostringstream of, of1;
	of<<parser.get(0)<<i<<"."<<"corr";
	out.open(of.str().c_str(), ios_base::out);
	if(!out) {
		cerr<<"Error opening output file "<<of.str()<<endl;
		exit(EXIT_FAILURE);
	};	
	
	//correct reads in this pool
	pTimer=new Timer("correct reads in pool");
	correctReads(i,numSearch,numHit);
	delete pTimer;
	//empty cache
	pools2bacs->makeEmpty();
	
	cout << "Pool " << i << ", numSearch " << numSearch << ", numHit " << numHit << endl;
	
	out.close();

}//end for i
	
	for (unsigned int r=0; r<MAX_POOL_SIZE; r++) {
		free(reads[r].header);
		free(reads[r].read);
	}
	free(reads);

return 0;

}

//}//namespace correction 
