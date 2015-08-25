#include <array>
#include <random>
#include <functional>
#include <unordered_map>

#include "Timer.h"
//#include "trie.h"
#include "classIO.h"

//#define DEBUG 1
#define MAX_DEPTH 10
#define MAX_KMERS(x)  (100 - x + 1)
#define MAX_POOL_SIZE 23000000

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
ofstream out, out2;
	
unsigned int numReads;
struct readFasta* reads;

struct equalTo {
	bool operator() (const set<unsigned short>& lhs, const set<unsigned short>& rhs) const {
	for (auto& b : lhs)
		if (rhs.find(b)==rhs.end())
			return false;
	return true;
	}
};

struct hashTo {
	size_t operator() (const set<unsigned short>& lhs) const {
	size_t hash=0;
	for (auto& b : lhs) {
		//if (b%2==0)
		//	hash=b<<1;
		hash += b;
	}	
	return (hash % 101);
	}
};

typedef pair<set<unsigned short>, vector<cset>> mypair;
typedef unordered_multimap<set<unsigned short>, vector<cset>, hashTo, equalTo> mymap;

void addCset(mymap& csets, cset& cst) {

if (cst.end - cst.begin + 1 < 3) 
	return;

cst.numPools=0;
short s = (cst.end - cst.begin + 1);
for (unsigned short i=0;i<POOLS;i++)
	if (cst.centroid[i]>0) {
		cst.numPools++;
		cst.centroid[i] /= s;
	};	

#if defined DEBUG
cout << "\t[addCset]centroid of newly added c-set between " << cst.begin << " and " << cst.end << endl;
cout << "\t";
for (unsigned short i=0;i<POOLS;i++)
	if (cst.centroid[i]>0) 
		printf("%d(%.2f) ", i + 1, cst.centroid[i]);
cout << endl;
#endif

if (csets.size()==0) {
	vector<cset> crt;
	crt.push_back(cst);
	//csets[cst.bacs]=crt;
	csets.insert(mypair(cst.bacs, crt));
	return;
}

bool found=false;
vector<cset> crt;
set<unsigned short> bacs, interBacs;
unordered_map<set<unsigned short>, vector<cset>, hashTo, equalTo>::iterator it1;
for (it1=csets.begin();it1!=csets.end();it1++) {
	bacs=it1->first; 
	crt=it1->second;
	
	if (cst.bacs.size()==0)
		interBacs=bacs;
	else {
		if (bacs.size()==0)
			interBacs=cst.bacs;
		else {
			interBacs.clear();
			for (auto& it2 : bacs) 
				if (cst.bacs.find(it2)!=cst.bacs.end()) 
					interBacs.insert(it2);
			
			if (interBacs.size()>0) {
				found=true;
				break;
			} else {	
				continue;
			}	
		};			
	};
	
	//find non-conflicting c-sets to join
	if (cst.begin>(short)(crt[crt.size() - 1].end + kmerSize)) {
		found=true;
		break;
	};	
}//end for

if (it1!=csets.end()) { 
	csets.erase(it1);
	crt.push_back(cst);
	//csets[interBacs]=crt;
	csets.insert(mypair(interBacs, crt));
};	

if (!found) {
	vector<cset> crt;
	crt.push_back(cst);
	//csets[cst.bacs]=crt;
	csets.insert(mypair(cst.bacs, crt));
};

};

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
	assert(crtKmer>=0);
	//make a working copy of the original k-mer
	//c=kmers[crtKmer];

}

double computeRho(double centroid[], unsigned short pools[], unsigned short flag) {

	double rho=0;
	double mean1=0, mean2=0;
	unsigned short numPools=0;
#if defined DEBUG
	if (flag)
		cout << "\t[computeRho]computing means over the following pools ";
#endif
	for (unsigned short p=0;p<POOLS;p++) {
		if (centroid[p]>0&&pools[p]>0) {
#if defined DEBUG
			if (flag)
				printf("%d(%.2f,%d) ", p + 1, centroid[p], pools[p]);
#endif
			mean1 += centroid[p];
			mean2 += pools[p];
			numPools++;
		};	
	};
#if defined DEBUG
if (flag)
	cout << endl;
#endif

	if (numPools<LOW2)
		return rho;
	
	mean1 /= numPools;
	mean2 /= numPools;

#if defined DEBUG
	if (flag)
		cout << "\t[computeRho]mean1 " << mean1 << ", mean2 " << mean2 << endl; 
#endif 

	double t1, t2;
	double s1=0, s2=0, s=0;
	for (unsigned int p=0;p<POOLS;p++) {
		if (centroid[p]>0&&pools[p]>0) {
			t1 = centroid[p] - mean1;
			t2 = pools[p] - mean2;
	
			s1 += pow(t1, 2);
			s2 += pow(t2, 2);
			s += t1 * t2;
		};	
	};

	if (s1==0) 
		s1=1;
	if (s2==0)
		s2=1;
	 
	 rho = s / (sqrt(s1)*sqrt(s2));

	 return rho;
};


bool correct(char* read, unsigned short direction, int crtKmer, int pos, int kmerPos, mybitsetx o, const vector<mybitsetx>& kmers, vector<mybitsetx>& kmersInHash, const vector<unsigned short>& orientation, vector<cset>& csets, unsigned short crtCset, unsigned int& numSearch, unsigned int& numHit) {

mysig sig;
mybitsetx g=kmers[crtKmer];
unsigned short orient=orientation[crtKmer];

if (g!=o) {
	mybitsetx c(o);
	c.invert();
	if (c<o)
		g=c;
	else
		g=o;
	//#pragma omp critical
	if (!ioh->searchCopy(g)) {
#if defined DEBUG
	cout << "[correct]kmer not found in the hashtable!" << endl;
#endif	
	return false;
	}
};

if (g.getNumPools()<LOW2) {
#if defined DEBUG
	cout << "[correct]too few pools " << g.getNumPools() << endl;
	for (unsigned short i=0;i<POOLS;i++) { 
		unsigned short gPoolCount = g.getPoolCount(i);
		if (gPoolCount>0) 
			cout << i+1 << "(" << gPoolCount << ") ";
	};
#endif	
	return false;
};

#if defined DEBUG
cout << endl;
cout << "[correct]Found, numPools " << g.getNumPools() << endl; 
for (unsigned short i=0;i<POOLS;i++) { 
	unsigned short gPoolCount = g.getPoolCount(i);
	if (gPoolCount>0) 
		cout << i+1 << "(" << gPoolCount << ") ";
};
cout << endl;	
cout << "[correct]crtCset " << crtCset << ", numPools " << csets[crtCset].numPools << ", centroid " << endl;
for (unsigned short i=0;i<POOLS;i++)  
	if (csets[crtCset].centroid[i]>0) 
		printf("%d(%.2f) ", i+1, csets[crtCset].centroid[i]);
cout << endl << endl;
#endif

//compute correlation coefficient
double rho=computeRho(csets[crtCset].centroid, g.getPools(), 1);

if (rho<=0.5)	{
#if defined DEBUG
	cout << "[correct]rho too small " << rho << endl;
#endif		
	return false;
};

#if defined DEBUG
cout << "[correct]rho " << rho << endl;
#endif		

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

};


void dfs(char* read, char* original, unsigned short direction, int crtKmer, unsigned short orient, int pos, int kmerPos, mybitsetx c, const vector<mybitsetx>& kmers, vector<mybitsetx>& kmersInHash, const vector<unsigned short>& orientation, vector<cset> csets, unsigned short crtCset, unsigned short depth, vector<unsigned short> posChanged, unsigned short& minCorrections, unsigned short& numCorrections,  unsigned int& numSearch, unsigned int& numHit) {
	

bool corrected=correct(read, direction, crtKmer, pos, kmerPos, c, kmers, kmersInHash, orientation, csets, crtCset,numSearch,numHit);

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
cout << endl << endl;	
//cout << "[dfs]kmer " << c << endl;
#endif	
		
#if defined DEBUG
cout << "[dfs]Trying " << c;
#endif	
	corrected=correct(read, direction, crtKmer, pos, kmerPos, c, kmers, kmersInHash, orientation, csets, crtCset, numSearch, numHit);
		
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


bool correctReads(unsigned int pool, unsigned int& numSearch, unsigned int& numHit) {

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

#pragma omp parallel default(none) shared(totalCorrect, noCset, canCorrect, cannotCorrect, totalChanged,totalUnchanged,cout,kmerSize,pool,ioh,trie,pools2bacs,reads,numReads,out,out2,numSearch,numHit,info) 
{

mysig sig;
bool gap=false;

ostringstream os;
ostringstream os2;

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
	unsigned short ePoolCount;
	mymap csets;
	
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
	ioh->searchCopy(e);
	
	kmers[k]=e;
	orientation[k]=orient;
	
	numBacs=0;
	numPools=e.getPools(pools);
	
	if (numPools>=LOW1 && numPools<=HIGH1) { //get bacs for this kmer 
		
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
		};
		
		//bacs=sig.getBacs();
		numBacs=sig.getNumBacs();	
		for (unsigned short i=0;i<numBacs;i++)
			bacs[i]=sig.getBac(i);
		
	}; //if (numPools>=LOW && numPools<=HIGH) 
	
#if defined DEBUG
	info[k].pos=k;
	info[k].numPools=numPools;
	info[k].numBacs=numBacs;
	for (unsigned int i=0;i<numBacs;i++)
		info[k].bacs[i]=bacs[i];
	for (unsigned int i=0;i<numPools;i++)
		info[k].pools[i]=pools[i];
#endif

	if (numPools<LOW1||(numPools>=LOW1&&numPools<=HIGH1&&numBacs==0)) {
		gap=true;
		continue;
	};	
	
	double rho=0;
	set<unsigned short> interBacs;
	if (crtCset.bacs.size()>0&&numBacs>0) {
	
		//should vote for same BACs
		for (auto& it : crtCset.bacs) {
			for (unsigned short b=0;b<numBacs;b++) 
				if (bacs[b]==it) {
					interBacs.insert(it);
					break;
				};
		};	
	} else if (crtCset.begin>=0) {
		
		//btw LOW2&&LOW1 or HIGH1&&HIGH2 compute rho
		rho=computeRho(crtCset.centroid, e.getPools(), 0);
#if defined DEBUG
		if (rho<=0.6)
			cout << "\t[correctReads]k " << k <<", rho too small " << rho << endl;
#endif		
	};
	
	if (k==(unsigned int)(crtCset.end + 1)&&(interBacs.size()>0||rho>0.6)) {//same cset	
		
		crtCset.end=k;
		if (interBacs.size()>0) 
			crtCset.bacs=interBacs;
		for (unsigned short p=0;p<POOLS;p++) {
			ePoolCount=e.getPoolCount(p);
			if (crtCset.centroid[p]>0&&ePoolCount>0)
				crtCset.centroid[p] += ePoolCount;
			else
				crtCset.centroid[p] = 0;
        	};
	} else {
		//add crtCset to  csets
		if (crtCset.begin>=0)
			addCset(csets, crtCset);
		
		//create new cset	
		crtCset.begin=k;
		crtCset.end=k;
		crtCset.bacs.clear();
		crtCset.bacs.insert(bacs, bacs + numBacs);	
		for (unsigned short p=0;p<POOLS;p++)
			crtCset.centroid[p]=e.getPoolCount(p);
#if defined DEBUG
		cout << "\t[correctReads]new cset @ k " << k << ", numPools " << numPools << endl;
#endif		
	}//end else
	
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
		cout << "(" << it2.begin << ", " << it2.end << ") ";
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
			os << 91 << endl; 
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
	unsigned short depth=0;
	while (depth<=MAX_DEPTH) {
#if defined DEBUG
cout << endl << "[correctReads]crt depth is " << depth << endl << endl;
#endif		
		correctable=false;
		for (it1=csets.begin();it1!=csets.end();it1++) {
			
			vector<cset> csts=it1->second;

			cst=0;
			findKmerToCorrect(read, csts, cst, orientation, direction, crtKmer, orient, pos, kmerPos);

			c=kmers[crtKmer];
#if defined DEBUG
cout << "[correctReads]crt csets: ";
for (auto& it2 : csts) 
	cout << "(" << it2.begin << ", " << it2.end << ") ";
cout << endl;	
cout << "[correctReads]Direction " << direction << " read pos to change " << pos << " orientation " << orient << " kmerPos to change " << kmerPos << endl;
cout << "[correctReads]read ";
for (int j=crtKmer;j<=(crtKmer+int(kmerSize)-1);j++)
	cout << read[j];
cout << endl << endl;	
//cout << "[correctReads]kmer " << c << endl;
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
				
				if (depth==0||(depth>0&&numCorrections==depth))
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
			os << 91 << endl; 
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
cout << "[correctReads]" << header << endl << endl;
#endif
		os << header << endl << original << endl;
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

#pragma omp critical (correct) 
out << os.str();

#pragma omp critical (uncorrected)
out2 << os2.str();

}//#pragma omp parallel 

cout << "[correctReads]noCset " << noCset << endl;
cout << "[correctReads]Total correct " << totalCorrect << endl;
cout << "[correctReads]Can correct " << canCorrect << endl;
cout << "[correctReads]Cannot correct " << cannotCorrect << endl;
cout << "[correctReads]Number of reads changed " << totalChanged << endl;
cout << "[correctReads]Number of reads unchanged " << totalUnchanged << endl;

if (cannotCorrect>0)
	return true;

return false;	

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

Parser parser;
char delimC[] = ".";
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

//unsigned short it;

ifstream in;
for (unsigned int i=0;i<91;i++) {
	
	//it=1;
	unsigned int numHit=0;
	unsigned int numSearch=0;
	
	ostringstream iss, oss, oss2; 
	iss<<parser.get(0)<<i<<"."<<parser.get(1);
	in.open(iss.str().c_str(), ios_base::in);
	if(!in) {
		cerr<<"Error opening file input file "<<iss.str()<<endl;
		exit(EXIT_FAILURE);
	};

	numReads = ioh->readReads(in, reads);
	printf("Number of reads read is %d\n", numReads);  
	
	//output files for corrected and uncorrected reads
	oss<<parser.get(0)<<i<<".corr";
	//oss2<<parser.get(0)<<i<<".uncorr"<<it;
	out.open(oss.str().c_str(), ios_base::out);
	//out2.open(oss2.str().c_str(), ios_base::out);
	
	if(!out) {
		cerr<<"Error opening output file "<<oss.str() << endl;
		exit(EXIT_FAILURE);
	};	
	//if(!out2) {
	//	cerr<<"Error opening output file "<<oss2.str() << endl;
	//	exit(EXIT_FAILURE);
	//};	
	
	pTimer=new Timer("correct reads in pool");
	//while (correctReads(i,numSearch,numHit)) { 
	//for (unsigned short n=0;n<1;n++) {	
		correctReads(i,numSearch,numHit);
		
		//in.close(); out2.close();
		
		//ostringstream iss2, oss3;
		//iss2<<parser.get(0)<<i<<".uncorr" << it;
		//oss3<<parser.get(0)<<i<< ".uncorr" << ++it;
		//in.open(iss2.str().c_str(), ios_base::in);
		//out2.open(oss3.str().c_str(), ios_base::out);
	
		//if(!in) {
		//	cerr<<"Error opening input file "<<iss2.str()<<endl;
		//	exit(EXIT_FAILURE);
		//};	
		//if(!out2) {
		//	cerr<<"Error opening output file "<<oss3.str()<<endl;
		//	exit(EXIT_FAILURE);
		//};	

		
		//numReads = ioh->readReads(in, reads);
		//printf("Number of reads read is %d\n", numReads);  
	//};
	delete pTimer;
	
	//empty cache
	//if (i%5==0)
	//	pools2bacs->makeEmpty();
	
	cout << "Pool " << i << ", numSearch " << numSearch << ", numHit " << numHit << endl;
	
	in.close();
	out.close();
	//out2.close();

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
