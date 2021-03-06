#include <random>
#include <array>
#include <functional>
#include <omp.h>

#include "Timer.h"
//#include "trie.h"
#include "classIO.h"

//#define DEBUG 1
#define MAX_DEPTH 4

//NEW IN THIS VERSION
//4. Multi-threading
//3. Estimating depth based on the distance between c-sets
//2. searchUpdate() only called at the end, after correcting the k-mers 
//1. When multiple c-sets are present which disagree with each other, treat them as correct one by one.

unsigned char c[8]={ 0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80};
unsigned char r[8]={ 0xFE, 0xFD, 0xFB, 0xF7, 0xEF, 0xDF, 0xBF, 0x7F};

unsigned int hashSizes[SIZES]={1000393, 2500009, 5000549, 10000643, 20000003, 30000001, 40000003, 50000017, 60000011,70000027, 80000023, 90000049, 100000567, 150000001, 200000477, 250000013, 300000641, 350000041, 400000651, 450000007, 500000461, 550000001, 600000607, 650000011, 700000747, 750000007, 800000533, 850000013, 900000637, 950000017, 1000000531, 1250000027, 1500000713, 1750000027, 2000000579, 2147483647}; 

using namespace std;
using namespace general;
using namespace IOFn;

namespace correction {

class MyCorrection {

private: 
	
	//all pools data
	IOF* ioh;
	unsigned int kmerSize;
	
	Trie* trie;
	HashTable<mysig>* pools2bacs;
	
	//per pool data
	string readFname;
	ofstream out;
	
	unsigned int pool;
	struct readFasta* reads;
	unsigned int numReads;
	
	//per read data
	vector<mybitsetx> kmers;
	vector<mybitsetx> kmersInHash;
	vector<unsigned short> orientation;
	
	unsigned int cores;

	unsigned int depth;
	unsigned int minCorrections;
	map<unsigned int, string> corrections;

	//per pool statistics
	unsigned int numSearch, numHit;

public:

MyCorrection(IOF* ioh, unsigned int kmerSize, Trie* trie, HashTable<mysig>* pools2bacs, string readFname, unsigned int pool, struct readFasta* reads, unsigned int numReads, unsigned int cores) {
	
	this->ioh=ioh;
	this->kmerSize=kmerSize;
	this->trie=trie;
	this->pools2bacs=pools2bacs;

	this->readFname=readFname;
	
	this->pool=pool;
	this->reads=reads;
	this->numReads=numReads;

	this->cores=cores;

	numSearch=0;
	numHit=0;
	
	//output files for corrected and uncorrected reads
	ostringstream of, of1;

	char delimC[] = ".";
	Parser parser;
	parser.update(delimC, readFname);

	of<<parser.get(0)<<pool<<"."<<"corr";
	out.open(of.str().c_str(), ios_base::out);
	if(!out) {
		cerr<<"Error opening output file "<<of.str()<<endl;
		exit(EXIT_FAILURE);
	};	

};

~MyCorrection() {
	cout << "Pool " << pool << ", numSearch " << numSearch << ", numHit " << numHit << endl;
	out.close();
}	


bool correct(char* read, int direction, int pos, int kmerPos, unsigned int crtCset, int crtKmer, mybitsetx c, unsigned short orient, vector<cset>& csets) {

bool correct=false;
//mybitsetx e=kmers[crtKmer];
//mybitsetx e=kmersInHash[crtKmer];
unsigned short bacs[MAX_BACS];
unsigned short pools[POOLS];
unsigned short numBacs, numPools; 

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
}	

numPools=g.getPools(pools);
if (numPools<LOW) {
#if defined DEBUG
	cout << "[correct]no pools found " << numPools << endl;
#endif	
	return false;	
}

if (numPools>HIGH)
	correct=true;

numBacs=0;
if (numPools>=LOW&&numPools<=HIGH) {
	mysig sig(pools, numPools);
	if (!pools2bacs->searchCopy(sig)) { 
		trie->search(pools, numPools, bacs, numBacs);
		
		numSearch++;
		
		sig.setBacs(bacs, numBacs);
		#pragma omp critical
		{
		pools2bacs->insert(sig);
		}
	} else {
		numHit++;
	}
	numBacs=sig.getNumBacs();	
	for (unsigned int i=0;i<numBacs;i++)
		bacs[i]=sig.getBac(i);
	
	if (numBacs==0)	{
#if defined DEBUG
		cout << "[correct]No bacs found " << numBacs << endl;
		cout << "[correct]numPools " << numPools << ": ";
		for (unsigned int i=0;i<numPools;i++)
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
		for (unsigned int i=0;i<numBacs;i++) {//does this kmer belong to crt cset?
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
		for (unsigned int i=0;i<numBacs;i++) 
			cout << bacs[i] << ",";
		cout << endl;
		return false;
#endif		
	}

	if (found) {
		set<unsigned short> newBacs;
		set<unsigned short> crtBacs=csets[crtCset].bacs;
		for (unsigned int i=0;i<numBacs;i++)
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
for (unsigned int i=0;i<numPools;i++) 
	cout << pools[i] << ",";
cout << endl;
cout << "[correct]numBacs " << numBacs << ": ";
for (unsigned int i=0;i<numBacs;i++) 
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
	
	//update pool counts in hashtable
	//ioh->searchUpdatePool(e, pool, -1);
	//ioh->searchUpdatePool(g, pool, 1);
	
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


void correctReads() {

#if defined DEBUG
struct kmerInfo info[MAX_KMERS];
#endif

unsigned short orient;
int crtKmer, pos, kmerPos; 

unsigned short pools[POOLS];
unsigned short numBacs, numPools;
unsigned short bacs[MAX_BACS];

omp_set_num_threads(cores);
printf("Number of cores is %d\n", cores);	

unsigned int r;
//unsigned int id;
//#pragma omp parallel private(id,r,kmers,kmersInHash,orientation,depth,minCorrections,corrections,numSearch,numHit) shared(pool,ioh,trie,pools2bacs,reads,numReads,out)
//{
//id = omp_get_thread_num();
//#pragma omp for schedule(dynamic)
for (r=0;r<numReads;r++) { 
	
	char* header = reads[r].header;
	char* read = reads[r].read;
	unsigned int num = strlen(read);
	char* original=(char*)malloc((num+1)*sizeof(char));
	strcpy(original, read);

	//count the number of errors the read contains
	char* errs=strchr(header, ':');
	unsigned int numErrors=0;
	while(*errs) if (*errs++ == '_') numErrors++;
	
	kmers.clear();
	orientation.clear();
	corrections.clear();
	
	unsigned int crtCset;
	vector<cset> csets;
	vector<cset> _csets;
	vector<unsigned short> posChanged;
	
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
			#pragma omp critical
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
			//crt.disagree=false;
			_csets.push_back(crt);	
#if defined DEBUG
			cout << "\t[new cset] at k " << k << ", numPools " << numPools << endl;
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
	} else {//create new cset	
		cset crt;
		crt.begin=k;
		crt.end=k;
		//crt.disagree=false;
		crt.bacs.insert(bacs, bacs + numBacs);	
		_csets.push_back(crt);	
#if defined DEBUG
		cout << "\t[new cset]k " << k << endl;
#endif		
	}	
	
	if (!found&&hasBacs) {
		conflicts=true;
		//for (auto it=_csets.begin();it!=_csets.end();it++) 
		//	if (it->bacs.size()>0)  
		//		it->disagree=true;
	}			

	}//end for k

#if defined DEBUG
cout << endl << header << " numErrors "<< numErrors << endl;
cout<< read << endl;
for (int k=0;k<MAX_KMERS;k++) {
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
next++;

depth=1;
while (crt!=_csets.end()) {
//for (auto it=_csets.begin();it!=_csets.end();it++) {

	if ((next!=_csets.end())&&(crt->end + (int)kmerSize + 1)<next->begin)
		depth++;
	
	if (conflicts&&crt->bacs.size()>0) {
		crt->end=crt->begin - 1;
		crt->bacs.clear();
	}	
	
	csets.push_back(*crt);

	prev=crt;
	crt=next;
	if (crt==_csets.end())
		break;
	next++;
}	

kmersInHash=kmers;

//use csets to fix errors by moving left and right

#if defined DEBUG
cout << "CORRECTION" << endl;
#endif
	crtCset=0;
	unsigned int oldCset=crtCset;
	unsigned int nextCset=crtCset + 1;

	int direction=0;
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
			if (nextCset!=csets.size())
				nextCset++;
			continue;
		}

		//update direction
		if ((csets[crtCset].begin - 1)!=prevECset) 
			direction=0;
		else if ((csets[crtCset].end + 1)!=nextBCset) 
			direction=1;
		break;
	}//end while
				
	//for correct reads don't do anything
	if (crtCset==csets.size()) { 
#if defined DEBUG
		cout << header << ":" << 1 << endl << endl;
#endif
		out << header << ":" << 1 << endl << read << endl;
		continue;
	};	
		
#if defined DEBUG
cout << endl << "[CORRECTION]crt cset begin " << csets[crtCset].begin << ", end " << csets[crtCset].end << " bacs ";
for (auto it=csets[crtCset].bacs.begin();it!=csets[crtCset].bacs.end();it++)
	cout << (*it) << ",";
cout << endl << endl; 
#endif	
	
	//determine pos in read&kmer which need to be changed 
	if (direction==0) {
		crtKmer=csets[crtCset].begin - 1;
		if (crtKmer==prevECset) {
			direction=1;
		} else {	
			pos=crtKmer;
			orient=orientation[crtKmer];
			if (orient)
				kmerPos=0;
			else
				kmerPos=kmerSize - 1;
		}		
	}
	if (direction==1) { 
		crtKmer=csets[crtCset].end + 1;
		pos=crtKmer + kmerSize - 1;
		orient=orientation[crtKmer];
		if (orient)
			kmerPos=kmerSize - 1;
		else
			kmerPos=0;
	}	

	e=kmers[crtKmer];
	
	//make a working copy of the original k-mer
	mybitsetx c(e);
	
#if defined DEBUG
cout << "crtKmer " << crtKmer << endl;
cout << "Direction " << direction << " read pos to change " << pos << " orientation " << orient << " kmerPos to change " << kmerPos << endl;

cout << "read ";
for (int j=crtKmer;j<=(crtKmer+int(kmerSize)-1);j++)
	cout << read[j];
cout << endl;	
cout << "kmer " << c;
#endif	

	//try all alternatives for first/last base of this k-mer
	vector<char> toTry;
	toTry.push_back(c.getBase(kmerPos));
	c.flip(2*kmerPos);
	toTry.push_back(c.getBase(kmerPos));
	c.flip(2*kmerPos);
	c.flip(2*kmerPos + 1);
	toTry.push_back(c.getBase(kmerPos));
	c.flip(2*kmerPos);
	toTry.push_back(c.getBase(kmerPos));
	
	minCorrections=100;
	while (depth<=MAX_DEPTH) {
#if defined DEBUG
cout << endl << "[correctReads]crt depth is " << depth << endl << endl;
#endif	
		for (unsigned int i=0;i<toTry.size();i++) { 
			c.setBase(kmerPos, toTry[i]);

#if defined DEBUG
cout << "[correctReads]crtKmer " << crtKmer << endl;
cout << "[correctReads]Trying " << c;
#endif	

			dfs(header, read, original, direction, pos, kmerPos, crtCset, prevECset, nextBCset, crtKmer, c, orient, csets, posChanged);
			
			//FIXME
			//Alternatively we can retrieve all solutions at the same depth 
			//and then no need to break out of the for loop here 
			if (minCorrections==depth)
				break;
	
		}//end for 
		
		if (minCorrections==depth)
			break;
		
		depth++;	
	}//end while
	
#if defined DEBUG
cout << "[correctReads]Num of solution paths " << corrections.size() << endl;
#endif

	if (corrections.size()==0) {

#if defined DEBUG
cout << "[correctReads]" << header << ":" << 0 << endl << endl;
#endif
		#pragma omp critical
		{
		out << header << ":" << 0 << endl << original << endl;
		}	
		continue;
	}//end if
	
#if defined DEBUG
cout << "[correctReads]Cost of best path is  " << minCorrections << endl;
#endif

#if defined DEBUG
cout << "[correctReads]" << header << ":" << 1 << endl << corrections[minCorrections] << endl << endl;
#endif
		#pragma omp critical
		{
		out << header << ":" << 1 << endl << corrections[minCorrections] << endl;
		}	

//update k-mer counts in the hashtable
for (unsigned int k=0; k<num-kmerSize+1; k++) {
	e=kmers[k];
	g=kmersInHash[k];
	if (g!=e) {
		#pragma omp critical
		{
		ioh->searchUpdatePool(e, pool, -1);
		ioh->searchUpdatePool(g, pool, 1);
		}
	}	
}

}//end for r
//}//#pragma omp parallel 
};

		
void dfs(char* header, char* read, char* original, int direction, int pos, int kmerPos, unsigned int crtCset, int prevECset, int nextBCset, int crtKmer, mybitsetx c, unsigned short orient, vector<cset> csets, vector<unsigned short> posChanged) {
	
bool corrected=correct(read, direction, pos, kmerPos, crtCset, crtKmer, c, orient, csets);
	
if (!corrected)
	return;
		
//record read pos which was succesfully changed
if (read[pos]!=original[pos])
	posChanged.push_back(pos);
	
unsigned int num = strlen(read);
unsigned int numCorrections=posChanged.size();
	
unsigned int oldCset=crtCset;
unsigned int nextCset=crtCset;
if (nextCset!=csets.size())
	nextCset++;
	
while (true) {
	
	oldCset=crtCset;
	
	//check stopping condition
	if (((csets[crtCset].begin - 1)==prevECset)&&((csets[crtCset].end + 1)==nextBCset)) { 
		crtCset=nextCset;
		if (nextCset!=csets.size())
			nextCset++;
	}

	while (crtCset!=csets.size()) {
		
		if (crtCset==oldCset)
			break;
			
		//update crt cset ends & neighbors	
		prevECset=csets[oldCset].end;
			
		nextBCset=num - kmerSize + 1;
		if (nextCset!=csets.size()) 
			nextBCset=csets[nextCset].begin;
				
		if (((csets[crtCset].begin - 1)==prevECset)&&((csets[crtCset].end + 1)==nextBCset)) {
			oldCset=crtCset;
			crtCset=nextCset;
			if (nextCset!=csets.size()) 
				nextCset++;
			continue;
		}
				
		//update direction
		if ((csets[crtCset].begin - 1)!=prevECset) 
			direction=0;
		else if ((csets[crtCset].end + 1)!=nextBCset) 
			direction=1;
		break;
		
	}//end while
	
	if (crtCset==csets.size())
		break;

	mybitsetx e;
	//determine pos in read&kmer which need to be changed next
	if (direction==0) {
		crtKmer=csets[crtCset].begin - 1;
		if (crtKmer==prevECset) {
#if defined DEBUG
cout << "\t[dfs]changing direction at " << (crtKmer + 1) << endl << endl;
#endif
			direction=1;
		} else {	
			pos=crtKmer;
			orient=orientation[crtKmer];
			if (orient)
				kmerPos=0;
			else
				kmerPos=kmerSize - 1;
		}		
	} 
	if (direction==1) { 
		crtKmer=csets[crtCset].end + 1;
		pos=crtKmer + kmerSize - 1;
		orient=orientation[crtKmer];
		if (orient)
			kmerPos=kmerSize - 1;
		else
			kmerPos=0;
	}	
		
	e=kmers[crtKmer];
		
	//make a copy of the original k-mer
	c=e;
	
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
		
	//FIXME
	//if (crtCset->disagree)
	//	break;

#if defined DEBUG
cout << "[dfs]Trying same k-mer " << c;
#endif	
	corrected=correct(read, direction, pos, kmerPos, crtCset, crtKmer, c, orient, csets);
		
	if (!corrected)
		break;
		
}//end while (true)		
	

if (crtCset==csets.size()) {
	if (numCorrections==depth) {

#if defined DEBUG
cout<<"[dfs]Adding NEW SOLUTION PATH with cost " << numCorrections <<endl<<endl;
#endif

		minCorrections=numCorrections;
		corrections[minCorrections]=read;
	} else {	
#if defined DEBUG
cout << "[dfs]numCorrections " << numCorrections << " RETURNING." << endl;
#endif
	}

	return;
}//end if
	
if ((numCorrections + 1)>depth) {
#if defined DEBUG
//cout << "[dfs]crtKmer " << crtKmer << endl;
cout<<"[dfs]Exceded MAX DEPTH. Returning."<<endl<<endl;
#endif
	return;
}

#if defined DEBUG
//cout<< endl << "[CORRECTION]crt cset begin " << csets[crtCset].begin << " end " << csets[crtCset].end << " bacs ";
//for (auto it=csets[crtCset].bacs.begin();it!=csets[crtCset].bacs.end();it++)
//	cout << (*it) << ",";
//cout << endl << endl; 
#endif	

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
	dfs(header, read, original, direction, pos, kmerPos, crtCset, prevECset, nextBCset, crtKmer, c, orient, csets, posChanged);
	if (minCorrections==depth)
		break;
	
};
	
};//end method definition 

};//end class definition

};//end namespace definition

