#include <unordered_map>
#include <random>
#include <array>
#include <functional>

#include "Timer.h"
//#include "trie.h"
#include "classIO.h"

//#define DEBUG 1
#define MAX_DEPTH 20
#define MAX_KMERS(x)  (100 - x + 1)
#define MAX_POOL_SIZE 2000000

//NEW IN THIS VERSION
//5. Cold/warm cache.
//4. Multi-threading
//3. Estimating depth based on the distance between c-sets
//2. searchUpdate() only called at the end, after correcting the k-mers 
//1. When multiple c-sets are present which disagree with each other, treat them as correct one by one.

unsigned char c[8]={ 0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80};
unsigned char r[8]={ 0xFE, 0xFD, 0xFB, 0xF7, 0xEF, 0xDF, 0xBF, 0x7F};

unsigned int hashSizes[SIZES]={5000549, 10000643, 20000003, 30000001, 40000003, 50000017, 60000011,70000027, 80000023, 90000049, 100000567, 150000001, 200000477, 250000013, 300000641, 350000041, 400000651, 450000007, 500000461, 550000001, 600000607, 650000011, 700000747, 750000007, 800000533, 850000013, 900000637, 950000017, 1000000531, 1250000027, 1500000713, 1750000027, 2000000579, 2147483647}; 

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

struct equalTo {
	bool operator() (const set<unsigned short>& lhs, const set<unsigned short>& rhs) const {
	for (auto it=lhs.begin();it!=lhs.end();it++)
		if (rhs.find(*it)==rhs.end())
			return false;
	return true;
	}
};

struct hashTo {
	size_t operator() (const set<unsigned short>& lhs) const {
	size_t hash=0;
	for (auto it=lhs.begin();it!=lhs.end();it++) {
		if ((*it)%2==0)
			hash=hash<<1;
		hash+=(*it);
	}	
	return hash;
	}
};

void addCset(unordered_map<set<unsigned short>, vector<cset>, hashTo, equalTo>& csets, cset& crtCset) {

if (crtCset.end - crtCset.begin + 1 < 3) 
	return;

if (csets.size()==0) {
	vector<cset> newCsets;
	newCsets.push_back(crtCset);
	csets[crtCset.bacs]=newCsets;
	return;
}

bool found=false;
set<unsigned short> newBacs;
unordered_map<set<unsigned short>, vector<cset>, hashTo, equalTo>::iterator it1;
for (it1=csets.begin();it1!=csets.end();it1++) {
	if (crtCset.bacs.size()==0) {
		found=true;
		it1->second.push_back(crtCset);

	} else {
		for (auto it2=it1->first.begin();it2!=it1->first.end();it2++) 
			if (crtCset.bacs.find(*it2)!=crtCset.bacs.end()) 
					newBacs.insert(*it2);
		if (it1->first.size()==0)
			newBacs=crtCset.bacs;
		if (newBacs.size()>0) {
			found=true;
			break; 
		}	
	}
}//end for

if (it1!=csets.end()) { 
	vector<cset> newCsets=it1->second;
	newCsets.push_back(crtCset);
	csets.erase(it1);
	csets[newBacs]=newCsets;
}	

if (!found) {
	vector<cset> newCsets;
	newCsets.push_back(crtCset);
	csets[crtCset.bacs]=newCsets;
}

}

void findKmerToCorrect(char* read, const vector<cset>& csets, unsigned short& crtCset, const vector<unsigned short>& orientation, unsigned short& direction, int& crtKmer, unsigned short& orient, int& pos, int& kmerPos) {

	unsigned int num=strlen(read);
	short oldCset=crtCset - 1;
	short nextCset=crtCset + 1;
	
	short prevECset=-1;
	if (oldCset>=0) 
		prevECset=csets[oldCset].end;
	
	short nextBCset=num - kmerSize + 1;
	if (nextCset<(short)csets.size()) 
		nextBCset=csets[nextCset].begin;
	
	if ((csets[crtCset].begin - 1)!=prevECset) 
		direction=0;
	else if ((csets[crtCset].end + 1)!=nextBCset) 
		direction=1;
	else 
		while (((csets[crtCset].begin - 1)==prevECset)&&((csets[crtCset].end + 1)==nextBCset)) {
			oldCset=crtCset;
			crtCset=nextCset;
			if (crtCset==(short)csets.size())
				break;
			nextCset++;
			
			prevECset=csets[oldCset].end;
			
			nextBCset=num - kmerSize + 1;
			if (nextCset!=(short)csets.size()) 
				nextBCset=csets[nextCset].begin;
			
			//update direction
			if ((csets[crtCset].begin - 1)!=prevECset) 
				{direction=0;break;}
			else if ((csets[crtCset].end + 1)!=nextBCset) 
				{direction=1;break;}

		}//end while
				
	if (crtCset==(short)csets.size()) 
		return;

#if defined DEBUG
//cout << endl << "[findKmerToCorrect]crt cset begin " << csets[crtCset].begin << ", end " << csets[crtCset].end << " bacs ";
//for (auto it=csets[crtCset].bacs.begin();it!=csets[crtCset].bacs.end();it++)
//	cout << (*it) << ",";
//cout << endl << endl; 
#endif	
	
	//determine pos in read&kmer which need to be changed 
	if (direction==0) {
		crtKmer=csets[crtCset].begin - 1;
		pos=crtKmer;
		orient=orientation[crtKmer];
		if (orient)
			kmerPos=0;
		else
			kmerPos=kmerSize - 1;
	} else { 
		crtKmer=csets[crtCset].end + 1;
		pos=crtKmer + kmerSize - 1;
		orient=orientation[crtKmer];
		if (orient)
			kmerPos=kmerSize - 1;
		else
			kmerPos=0;
	}	

	//make a working copy of the original k-mer
	//c=kmers[crtKmer];

}

bool correct(char* read, unsigned short direction, int crtKmer, int pos, int kmerPos, mybitsetx c, vector<mybitsetx>& kmersInHash, const vector<unsigned short>& orientation, vector<cset>& csets, unsigned short crtCset, unsigned int& numSearch, unsigned int& numHit) {

mysig sig;
bool correct=false;
unsigned short pools[POOLS];
unsigned short bacs[MAX_BACS];
unsigned short numBacs=0, numPools=0; 
unsigned short orient=orientation[crtKmer];

mybitsetx g;
mybitsetx o(c);
c.invert();
if (c<o)
	g=c;
else
	g=o;

//#pragma omp critical
//{
if (!ioh->searchCopy(g)) {
//}
#if defined DEBUG
	cout << "[correct]kmer not found in the hashtable!" << endl;
#endif	
	return false;
};	

numBacs=0;
numPools=g.getPools(pools);
//pools=g.getPools();
//numPools=g.getNumPools();

if (numPools<LOW) {
#if defined DEBUG
	cout << "[correct]No pools found " << numPools << endl;
#endif	
	return false;	
};

if (numPools>HIGH)
	correct=true;

if (numPools>=LOW&&numPools<=HIGH) {
	
	sig.setPools(pools, numPools);
	
	if (pools2bacs->searchCopy(sig)) { 

		#pragma omp atomic
		numHit++;

	} else {	
		sig.resetNumBacs();
		trie->search(pools, numPools, sig.getBacs(), sig.getNumBacs());
		
		#pragma omp critical (insert)
		{
		pools2bacs->insert(sig);
		}
		
		#pragma omp atomic
		numSearch++;
	} 
	//bacs=sig.getBacs();
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
	};
	
	bool hasBacs=false, found=false;
	//check that bacs found agree with other cset bacs
	//for (auto it=csets.begin();it<csets.begin() + crtCset + 1;it++) {
	for (auto it=csets.begin() + crtCset;it>=csets.begin();it--) {
		if (it->bacs.size()==0)
			continue;
			
		hasBacs=true;
		
		for (auto& it2 : it->bacs) {
			for (unsigned short i=0;i<numBacs;i++)
				if (it2==bacs[i]) {
					found=true;
					break;
				}
			if (found)
				break;
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
		for (auto& it : crtBacs) 
			for (unsigned short i=0;i<numBacs;i++)
				if (it==bacs[i]) {
					newBacs.insert(it);
					break;
				}		
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

void dfs(char* read, char* original, unsigned short direction, int crtKmer, unsigned short orient, int pos, int kmerPos, mybitsetx c, const vector<mybitsetx>& kmers, vector<mybitsetx>& kmersInHash, const vector<unsigned short>& orientation, vector<cset> csets, unsigned short crtCset, unsigned short depth, vector<unsigned short> posChanged, unsigned short& minCorrections, unsigned short& numCorrections,  unsigned int& numSearch, unsigned int& numHit) {
	

bool corrected=correct(read, direction, crtKmer, pos, kmerPos, c, kmersInHash, orientation, csets, crtCset,numSearch,numHit);

if (!corrected)
	return;

//record read pos which was succesfully changed
if (read[pos]!=original[pos])
	posChanged.push_back(pos);

numCorrections=posChanged.size();

while (true) {
	
	findKmerToCorrect(read, csets, crtCset, orientation, direction, crtKmer, orient, pos, kmerPos);
		
	if (crtCset==csets.size())
		break;

	//make a copy of the original k-mer
	c=kmers[crtKmer];
	
#if defined DEBUG
cout << "[dfs]crtKmer " << crtKmer << endl;
cout << "[dfs]num corrections so far " << numCorrections << endl;
if (posChanged.size()>0)
	cout << "[dfs]Updating kmer with existing read changes" << endl;
#endif

	int kp;
#if defined DEBUG		
	char oldBase;
#endif			
	//update crt kmer with possible base changes from previous fixed kmers
	for (auto& p : posChanged) {
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
cout << "[dfs]crtKmer " << crtKmer << endl;
cout << "[dfs]Direction " << direction << " read pos to change " << pos << " orientation " << orient << " kmerPos to change " << kmerPos << endl;
cout << "[dfs]read ";
for (int j=crtKmer;j<=(crtKmer+int(kmerSize)-1);j++)
	cout << read[j];
cout << endl;	
cout << "[dfs]kmer " << c;
#endif	
		
#if defined DEBUG
cout << "[dfs]Trying k-mer " << c;
#endif	
	corrected=correct(read, direction, crtKmer, pos, kmerPos, c, kmersInHash, orientation, csets, crtCset, numSearch, numHit);
		
	if (!corrected)
		break;
		
}//end while (true)		
	

if (crtCset==csets.size()) {
#if defined DEBUG
cout<<"[dfs]Found SOLUTION PATH with cost " << numCorrections <<endl<<endl;
#endif
	minCorrections=numCorrections;
	return;
}//end if
	
if ((numCorrections + 1)>depth) {
#if defined DEBUG
cout<<"[dfs]Exceded MAX DEPTH. Returning."<<endl<<endl;
#endif
	return;
}

//recurse
array<char,3> toTry;
c.flip(2*kmerPos);
toTry[0]=c.getBase(kmerPos);
c.flip(2*kmerPos);
c.flip(2*kmerPos + 1);
toTry[1]=c.getBase(kmerPos);
c.flip(2*kmerPos);
toTry[2]=c.getBase(kmerPos);

//bool correctable =false;
//unsigned short oldNumCorrections=numCorrections;
for (unsigned int i=0;i<toTry.size();i++) { 
	c.setBase(kmerPos, toTry[i]);
#if defined DEBUG
cout << "[dfs]crtKmer " << crtKmer << endl;
cout << "[dfs]Trying " << c;
#endif	
	dfs(read, original, direction, crtKmer, orient, pos, kmerPos, c, kmers, kmersInHash, orientation, csets, crtCset, depth, posChanged, minCorrections, numCorrections, numSearch, numHit);
	
	//if (numCorrections>oldNumCorrections)
	//	correctable;
	if (minCorrections==depth)
		break;
};

//if (correctable)
//	return true;

//return false;	

};//end method definition 


void correctReads(unsigned int pool, unsigned int& numSearch, unsigned int& numHit) {

//#if defined DEBUG
struct kmerInfo info[MAX_KMERS(kmerSize)];
//#endif

omp_set_num_threads(cores);
printf("Number of cores is %d\n", cores);	

unsigned int totalCorrect=0;
unsigned int noCset=0;
unsigned int canCorrect=0;
unsigned int cannotCorrect=0;
unsigned int totalChanged=0;
unsigned int totalUnchanged=0;

#pragma omp parallel default(none) shared(totalCorrect, noCset, canCorrect, cannotCorrect, totalChanged,totalUnchanged,cout,kmerSize,pool,ioh,trie,pools2bacs,reads,numReads,out,numSearch,numHit,info) 
{

mysig sig;
bool gap=false;

ostringstream os;

unsigned short orient;
int crtKmer=0, pos=0, kmerPos=0; 

unsigned short numPools=0;
unsigned short pools[POOLS];
unsigned short numBacs=0;
unsigned short bacs[MAX_BACS];

vector<mybitsetx> kmers(MAX_KMERS(kmerSize));
vector<mybitsetx> kmersInHash(MAX_KMERS(kmerSize));
vector<unsigned short> orientation(MAX_KMERS(kmerSize));

#pragma omp for schedule(static)
for (unsigned int r=0;r<numReads;r++) { 
	
	char* header = reads[r].header;
	//out << header << endl;
	char* read = reads[r].read;
	unsigned int num = strlen(read);
	char* original=(char*)malloc((num+1)*sizeof(char));
	strcpy(original, read);

	//char* errs=strchr(header, ':');
	//unsigned int numErrors=0;
	//while(*errs) if (*errs++ == '_') numErrors++;
	
	gap=false;
	cset crtCset;
	unordered_map<set<unsigned short>, vector<cset>,hashTo, equalTo> csets;

	class mybitsetx b,c,d,e,g;
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
	kmers[k]=e;
	orientation[k]=orient;
	
	numBacs=0;
	numPools=e.getPools(pools);
	//pools=e.getPools();
	//numPools=e.getNumPools();
	
	if (numPools>=LOW && numPools<=HIGH) { //get bacs for this kmer 
		
		sig.setPools(pools, numPools);
		
		if (pools2bacs->searchCopy(sig)) { 
			
			#pragma omp atomic
			numHit++;
		
		} else {	
			sig.resetNumBacs();
			trie->search(pools, numPools, sig.getBacs(), sig.getNumBacs());
			
			#pragma omp critical (insert)
			{
			pools2bacs->insert(sig);
			}
			
			#pragma omp atomic
			numSearch++;
		}
		
		//bacs=sig.getBacs();
		numBacs=sig.getNumBacs();	
		for (unsigned short i=0;i<numBacs;i++)
			bacs[i]=sig.getBac(i);
		
	} //if (numPools>=LOW && numPools<=HIGH) 
	
#if defined DEBUG
	info[k].pos=k;
	info[k].numPools=numPools;
	info[k].numBacs=numBacs;
	for (unsigned int i=0;i<numBacs;i++)
		info[k].bacs[i]=bacs[i];
	for (unsigned int i=0;i<numPools;i++)
		info[k].pools[i]=pools[i];
#endif

	if (numPools<LOW||(numPools>=LOW&&numPools<=HIGH&&numBacs==0)) {
		gap=true;
		continue;
	}	
	
	if (numPools>HIGH) { 
		if (crtCset.bacs.size()==0&&k==(unsigned int)(crtCset.end + 1)) { 
			crtCset.end=k;
		} else {
			//add crtCset to  csets
			if (crtCset.begin>=0)
				addCset(csets, crtCset);
			//create new cset
			crtCset.begin=k;
			crtCset.end=k;
			crtCset.bacs.clear();
#if defined DEBUG
		cout << "\t[correctReads]new cset @ k " << k << ", numPools " << numPools << endl;
#endif		
		};	
		
		continue;
	} //end if (numPools>HIGH)	
	
	set<unsigned short> newBacs;
	for (auto& it : crtCset.bacs)
		for (unsigned short index=0;index<numBacs;index++) 
			if (bacs[index]==it) {
				newBacs.insert(it);
				break;
			}	
	
	if (newBacs.size()>0&&k==(unsigned int)(crtCset.end + 1)) {//same cset	
		crtCset.end=k;
		crtCset.bacs=newBacs;
	} else {
		//add crtCset to  csets
		if (crtCset.begin>=0)
			addCset(csets, crtCset);
		//create new cset	
		crtCset.begin=k;
		crtCset.end=k;
		crtCset.bacs.clear();
		crtCset.bacs.insert(bacs, bacs + numBacs);	
#if defined DEBUG
		cout << "\t[correctReads]new cset @ k " << k << ", numPools " << numPools << endl;
#endif		
	};//end else
	
	}//end for k

//add last cset
if (crtCset.begin>=0)
	addCset(csets, crtCset);

#if defined DEBUG
cout << endl << header << ", numKmers " << (num - kmerSize + 1) << endl << read << endl;
for (unsigned int k=0;k<num-kmerSize+1;k++) {
	cout << "kmer " << info[k].pos << " numPools " << info[k].numPools << " numBacs " << info[k].numBacs << " bacs ";
	if (info[k].numBacs<=21) { 
		for (unsigned int i=0;i<info[k].numBacs;i++) 
			cout << info[k].bacs[i] << ",";
	}
	if (info[k].numPools>=3 && info[k].numPools<=14) {
		cout << " pools ";
		for (unsigned int i=0;i<info[k].numPools;i++) 
			cout << info[k].pools[i] << ",";
	}
	cout << endl;	
}
cout << header << endl << endl;
#endif

kmersInHash=kmers;

#if defined DEBUG
unsigned short i=0;
cout << "csets" << endl;
for (auto& it : csets) {
	cout << "cset set " << i++ << endl;
	cout << "bacs: ";
	for (auto& it2 : it.first)
		cout << it2 << " ";
	cout << endl << "csets: ";
	for (auto& it2 : it.second) 
		cout << "(" << it2.begin << ", " << it2.end << "), ";
	cout << endl;	
}
#endif

#if defined DEBUG
cout << "CORRECTION" << endl;
#endif
	
	if (csets.size()==0) {
		os << header << endl << read << endl;
		os << "-1 -1 -1" << endl;
		os << 0 << endl;
		
		#pragma omp atomic
		noCset++;

		continue;
	}

	if (!gap) {
#if defined DEBUG
		cout << "[correctReads]Read already correct" << endl;
		cout <<"[correctReads]" << header << endl << endl;
#endif
		set<unsigned short> bacs2=csets.begin()->first;
		os << header << endl << read << endl;
		if (bacs2.size()==0||bacs2.size()>SPARSITY) {
			os << "-1 -1 -1" << endl;
			os << 0 << endl; 
		} else {
			for (auto& bac : bacs2) 
				os << bac << " ";
			for (unsigned short j=bacs2.size();j<3;j++) 
				os << (-1) << " ";
			os << endl;
		}
		
		#pragma omp atomic
		totalCorrect++;

		continue;
	}
	
	unsigned short cst=0;
	unsigned short direction=0;
	
	array<char,4> toTry;
	unsigned short numCorrections=0;
	unsigned short minCorrections=100;
	vector<unsigned short> posChanged;
	
	unordered_map<set<unsigned short>,vector<cset>, hashTo, equalTo>::iterator it1;
	
	bool correctable=false;
	unsigned short depth=1;
	while (depth<=MAX_DEPTH) {
#if defined DEBUG
cout << endl << "[correctReads]crt depth is " << depth << endl << endl;
#endif		
		for (it1=csets.begin();it1!=csets.end();it1++) {
			
			vector<cset> csts=it1->second;

			cst=0;
			findKmerToCorrect(read, csts, cst, orientation, direction, crtKmer, orient, pos, kmerPos);

			c=kmers[crtKmer];
#if defined DEBUG
cout << "[correctReads]crtKmer " << crtKmer << endl;
cout << "[correctReads]Direction " << direction << " read pos to change " << pos << " orientation " << orient << " kmerPos to change " << kmerPos << endl;
cout << "[correctReads]read ";
for (int j=crtKmer;j<=(crtKmer+int(kmerSize)-1);j++)
	cout << read[j];
cout << endl;	
cout << "[correctReads]kmer " << c;
#endif	
		
			//try all alternatives for first/last base of this k-mer
			toTry[0]=c.getBase(kmerPos);
			c.flip(2*kmerPos);
			toTry[1]=c.getBase(kmerPos);
			c.flip(2*kmerPos);
			c.flip(2*kmerPos + 1);
			toTry[2]=c.getBase(kmerPos);
			c.flip(2*kmerPos);
			toTry[3]=c.getBase(kmerPos);
		
			for (unsigned int i=0;i<toTry.size();i++) { 
				c.setBase(kmerPos, toTry[i]);

#if defined DEBUG
cout << "[correctReads]crtKmer " << crtKmer << endl;
cout << "[correctReads]Trying " << c;
#endif	

				dfs(read, original, direction, crtKmer, orient, pos, kmerPos, c, kmers, kmersInHash, orientation, csts, cst, depth, posChanged, minCorrections, numCorrections, numSearch, numHit);
			
				if (minCorrections==depth)
					break;
				
				if (numCorrections>0)
					correctable=true;

			}//end for 
		
			if (minCorrections==depth)
				break;

		}//end for	
		
		if (minCorrections==depth)
			break;

		if (!correctable)
			break;
		
		depth++;	
	
	}//end while	
	
	if (minCorrections==depth) {

#if defined DEBUG
cout << "[correctReads]Cost of best path is  " << minCorrections << endl;
cout << "[correctReads]" << header << ":" << 1 << endl << read << endl << endl;
#endif
		os << header << endl << read << endl;
		if (it1->first.size()==0||it1->first.size()>SPARSITY) {
			os << "-1 -1 -1" << endl;
			os << 0 << endl; 
		} else {
			for (auto& bac : it1->first) 
				os << bac << " ";
			for (unsigned short j=it1->first.size();j<3;j++) 
				os << (-1) << " ";
			os << endl;
		}
		
		#pragma omp atomic
		canCorrect++;
	} else {

#if defined DEBUG
cout << "[correctReads]Unable to correct " << endl;
cout << "[correctReads]" << header << ":" << 0 << endl << endl;
#endif
		os << header << endl << read << endl;
		os << "-1 -1 -1" << endl;
		os << 0 << endl; 
		#pragma omp atomic
		cannotCorrect++;
		if (numCorrections>0) {
			#pragma omp atomic
			totalChanged++;
		} else {
			//DEBUG
			//os << header << endl << read << endl;
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

};

int main(int argc, char* argv[]) {

cout << "size of unsigned " << sizeof(unsigned) << endl;
cout << "size of unsigned long " << sizeof(unsigned long) << endl;
cout << "size of unsigned long long " << sizeof(unsigned long long) << endl;

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
pools2bacs->allocHash(2);

//reads=(struct readFasta*)malloc(MAX_POOL_SIZE*sizeof(struct readFasta));
reads=new struct readFasta[MAX_POOL_SIZE];
if (!reads) {
	perror("Error allocating the reads!");
	exit(1);
}
/*for (unsigned int r=0; r<MAX_POOL_SIZE; r++) {
	reads[r].header=NULL;
	reads[r].read=NULL;
}*/

ifstream in;
for (unsigned int i=0;i<91;i++) {
	
	unsigned int numHit=0;
	unsigned int numSearch=0;
	
	ostringstream iss, oss;
	iss<<parser.get(0)<<i<<"."<<parser.get(1);
	in.open(iss.str().c_str(), ios_base::in);
	if(!in) {
		cerr<<"Error opening intput file "<<iss.str()<<endl;
		exit(EXIT_FAILURE);
	};	
	
	numReads = ioh->readReads(in, reads);
	printf("Number of reads read is %d\n", numReads);  
	
	//output files for corrected and uncorrected reads
	oss<<parser.get(0)<<i<<".corr";
	out.open(oss.str().c_str(), ios_base::out);
	if(!out) {
		cerr<<"Error opening output file "<<oss.str()<<endl;
		exit(EXIT_FAILURE);
	};	
	
	//correct reads in this pool
	pTimer=new Timer("correct reads in pool");
	correctReads(i,numSearch,numHit);
	delete pTimer;
	//empty cache
	//if (i%5==0)
	//	pools2bacs->makeEmpty();
	
	cout << "Pool " << i << ", numSearch " << numSearch << ", numHit " << numHit << endl;
	
	in.close();
	out.close();

}//end for i
	
	/*for (unsigned int r=0; r<MAX_POOL_SIZE; r++) {
		free(reads[r].header);
		free(reads[r].read);
	}
	free(reads);*/
	delete [] reads;

return 0;

}

//}//namespace correction 
