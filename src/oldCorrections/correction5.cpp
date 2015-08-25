#include <omp.h>
#include <map>
#include <set>
#include <random>
#include <functional>

#include "Timer.h"
#include "trie.h"
#include "classIO.h"

//#define DEBUG 1

unsigned char c[8]={ 0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80};
unsigned char r[8]={ 0xFE, 0xFD, 0xFB, 0xF7, 0xEF, 0xDF, 0xBF, 0x7F};

unsigned int hashSizes[SIZES]={1000393, 2500009, 5000549, 10000643, 20000003, 30000001, 40000003, 50000017, 60000011,70000027, 80000023, 90000049, 100000567, 150000001, 200000477, 250000013, 300000641, 350000041, 400000651, 450000007, 500000461, 550000001, 600000607, 650000011, 700000747, 750000007, 800000533, 850000013, 900000637, 950000017, 1000000531, 1250000027, 1500000713, 1750000027, 2000000579, 2147483647}; 

using namespace std;
using namespace general;
using namespace IOFn;

Trie trie;
class IOF* p_ref;
//int numBacs;
//int numPools;
//unsigned short bacs[7];
//unsigned short pools[POOLS];

//typedef map<int, set<unsigned short>> pos2bacs;
//typedef map<int, mybitsetx> pos2kmers;
//typedef map<int, unsigned short> pos2orientation;

namespace IOFn {

//check intermediate k-mers between crtBrk and nextBrk -once kmer corresposnding to crtBrk was succesfully error-corrected-to make sure all k-mers are now valid and detect hidden brks 
//pos is the read pos which was succesfully changed
void IOF::checkIntermediateKmers(int pool, char* read, const vector<int>& posChanged, int direction, int begin, int end, struct cset* crt, struct cset* neighbor, const vector<mybitsetx>& kmers, const vector<unsigned short>& orientation, HashTable<mysig>& pools2bacs, unsigned int& numSearch, unsigned int& numHit) {

bool found=false;

unsigned short bacs[MAX_BACS];
unsigned short pools[POOLS];
unsigned int numBacs, numPools; 

//check all k-mers between crt brk and the next one
//because there might be clustered errors which obfuscated some brks
int p=begin;
while (true) {
	if (direction==0)
		p--;
	else 
		p++;
	if (p==end)
		break;
	
	mybitsetx e=kmers[p];
	unsigned short orient=orientation[p];
	
	//update pool counts in hashtable
	sh.searchUpdatePool(e, pool, -1);
	
#if defined DEBUG
	//if (posChanged.size()>0)
		//cout << "[checkIntermediateKmers]Updating kmer with existing read changes" << endl;
#endif
	int pos, kmerPos;
	//updating crtBrk kmer with possible base changes from previous fixed brk
	for (auto it=posChanged.begin();it!=posChanged.end();it++) {
		pos=(*it);
		if ((p<=pos)&&(pos<=(p+int(kmerSize)-1))) {
			//char oldBase;
			if (orient) {
				kmerPos=pos-p;
				//oldBase=e.getBase(kmerPos);
				e.setBase(kmerPos, read[pos]);
			} else {
				kmerPos=kmerSize-1-(pos-p);
				//oldBase=e.getBase(kmerPos);
				e.setBase(kmerPos, mybitset::rev(read[pos]));
			}	
#if defined DEBUG
			//cout<<"\t[checkIntermediateKmers]kmer p "<<p<<" read pos prev changed "<<pos<<" new read base "<<read[pos]<< " kmerPos changed "<<kmerPos<<" old kmer base "<<oldBase<<" new kmer base "<<e.getBase(kmerPos)<<" orientation "<<orient<<endl;   
#endif

	}//end if
	}//end for it	
	
	mybitsetx g;
	mybitsetx o(e);
	e.invert();
	if (o<e)
		g=o;
	else
		g=e;

	found=false;
	numBacs=0; numPools=0;
	if (sh.searchCopy(g)) {//search kmer in hashtable
		
		//update pool counts in hashtable
		sh.searchUpdatePool(g, pool, 1);
		
		numPools=g.count(pools);
		
		if (numPools>MAXREP) 
			found=true;
		
		if (numPools>=4&&numPools<=MAXREP) {
			mysig sig(pools, numPools);
			if (!pools2bacs.searchCopy(sig)) { 
				trie.search(pools, numPools, bacs, numBacs);
				numSearch++;
				sig.setBacs(bacs, numBacs);
				pools2bacs.insert(sig);
			} else {
				numHit++;
			}
			numBacs=sig.getNumBacs();	
			for (unsigned int i=0;i<numBacs;i++)
				bacs[i]=sig.getBac(i);
			
			for (unsigned int i=0;i<numBacs;i++) {//check crt cset
				if (crt->bacs.find(bacs[i])!=crt->bacs.end()) {
					found=true;
					break;
				}
			}
		
			if (!found) {//check next cset too
				if (neighbor!=NULL) {//neighbor can be either prev or next cset
					for (unsigned int i=0;i<numBacs;i++) {
						if (neighbor->bacs.find(bacs[i])!=neighbor->bacs.end()) {	
							found=true;
							break;
						}
					}
				}
			}	
				
		}
	}	
	
#if defined DEBUG
	cout << "\t[checkIntermediate]checking p " << p << " direction " << direction << " numPools " << numPools << " numBacs " << numBacs << endl;
#endif		

	if (!found) {//found new error
#if defined DEBUG
		cout << "\t[checkIntermediateKmers]checked kmer at pos  " << p << endl; 
		cout << "\t[checkIntermediateKmers]numPools " << numPools << endl;
		cout << "\t[checkIntermediateKmers]numBacs " << numBacs << ": ";
		for (unsigned int i=0;i<numBacs;i++) 
			cout << bacs[i] << ",";
		cout << endl;	
		cout << "\t[checkIntermediateKmers]kmer to correct next " << p << endl;
#endif		
		break;
	}
	if (direction==0)
		crt->begin--;
	else	
		crt->end++;
}//while (p!=end) {

};


bool IOF::correct(int pool, char* read, int direction, int pos, int kmerPos, struct cset* crt, struct cset* neighbor, const mybitsetx& e, mybitsetx c, unsigned short orient, HashTable<mysig>& pools2bacs, unsigned int& numSearch, unsigned int& numHit) {

unsigned int i;
bool correct=false;
unsigned short bacs[MAX_BACS];
unsigned short pools[POOLS];
unsigned int numBacs, numPools; 

mybitsetx g;
mybitsetx o(c);
c.invert();
if (c<o)
	g=c;
else
	g=o;

//DEBUG
//cout << "[correct]Searching " << g; 
if (!sh.searchCopy(g)) {
#if defined DEBUG
	cout << "kmer not found in the hashtable!" << endl;
#endif	
	return false;
}	

numPools=g.count(pools);
if (numPools<4) {
#if defined DEBUG
	cout << "[correct]no pools found " << numPools << endl;
#endif	
	return false;	
}

if (numPools>MAXREP)
	correct=true;

numBacs=0;
if (numPools>=4&&numPools<=MAXREP) {
	mysig sig(pools, numPools);
	if (!pools2bacs.searchCopy(sig)) { 
		trie.search(pools, numPools, bacs, numBacs);
		numSearch++;
		sig.setBacs(bacs, numBacs);
		pools2bacs.insert(sig);
	} else {
		numHit++;
	}
		numBacs=sig.getNumBacs();	
		for (unsigned int i=0;i<numBacs;i++)
			bacs[i]=sig.getBac(i);
	
	if (numBacs==0)	{
#if defined DEBUG
		cout << "[correct]no bacs found " << numBacs << endl;
#endif		
		return false;
	}
	
	//check that bacs found agree with crt cset bacs
	bool found=false;
	for (i=0;i<numBacs;i++) {
		if (crt->bacs.find(bacs[i])!=crt->bacs.end()) {	
			found=true;
			break;
		}
	}
	
	if (!found) {//check next cset too
		if (neighbor!=NULL) {//neighbor can be either prev or next cset
			for (i=0;i<numBacs;i++) {
				if (neighbor->bacs.find(bacs[i])!=neighbor->bacs.end()) {	
					found=true;
					break;
				}
			}
		}
	}	

	if (found) {
		correct=true;
		for (i=0;i<numBacs;i++) 
			crt->bacs.insert(bacs[i]);
	}	
	
}//if (numPools>=4&&numPools<=MAXREP) {

#if defined DEBUG
cout << "[correct]Found " << g; 
cout << "[correct]numPools " << numPools << endl;
cout << "[correct]numBacs " << numBacs << ": ";
for (i=0;i<numBacs;i++) 
	cout << bacs[i] << ",";
cout << endl;
#endif

if (correct) {
	//extend crt cset
	if (direction==0)
		crt->begin--;
	else	
		crt->end++;
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
	//update pool counts in hashtable
	sh.searchUpdatePool(e, pool, -1);
	sh.searchUpdatePool(g, pool, 1);
	
#if defined DEBUG
	cout << "[correct]read pos " << pos << " old base " << oldBase << " new base " << newBase << endl;
#endif
	return true;
}
	
return false;

};


void IOF::correctReads(int pool, struct readFasta* reads, int numReads, HashTable<mysig>& pools2bacs)
{

//output files for corrected and uncorrected reads
ostringstream of, of1;

char delimC[] = ".";
Parser parser;
parser.update(delimC, ReadFname);

ofstream out;
of<<parser.get(0)<<pool<<"."<<"corr";
out.open(of.str().c_str(), ios_base::out);
if(!out) {
	cerr<<"Error opening output file "<<of.str()<<endl;
	exit(EXIT_FAILURE);
}	

#if defined DEBUG
unsigned int i;
struct kmerInfo info[MAX_KMERS];
#endif

unsigned int numSearch=0;
unsigned int numHit=0;
for (int r=0;r<numReads;r++) { 
	
	char* header = reads[r].header;
	char* read = reads[r].read;
	int num = strlen(read);
	char* original=(char*)malloc((num+1)*sizeof(char));
	strcpy(original, read);

	//count the number of errors the read contains
	char* errs=strchr(header, ':');
	int numErrors=0;
	while(*errs) if (*errs++ == '_') numErrors++;
	
	vector<mybitsetx> kmers;
	vector<cset> csets;
	unsigned short crtCset;
	unsigned int numPools;
	unsigned short pools[POOLS];
	unsigned short orient;
	vector<unsigned short> orientation;

	class mybitsetx b,d,e,g;
	for (unsigned int k=0; k<num-kmerSize+1; k++) {
		int z=0;
		int j=kmerSize*2-1;
		//bool findN=false;
		unsigned int c;
		for (c=0; c<kmerSize; c++) {
			switch (read[c+k]) {
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
	
	bool inHash=sh.searchCopy(e); 
	
	kmers.push_back(e);
	orientation.push_back(orient);
	
	//get positive pools and their number for this k-mer
	numPools=e.count(pools);

	unsigned int numBacs=0;
	unsigned short bacs[MAX_BACS];
	//unsigned short bacs[3];
	//for (unsigned int j=0;j<7;j++)
	//	cout << bacs[j] << " ";
	//cout << endl;	
	if (numPools>=4 && numPools<=MAXREP) { //get bacs for this kmer 
		mysig sig(pools, numPools);
		if (!pools2bacs.searchCopy(sig)) { 
			trie.search(pools, numPools, bacs, numBacs);
			numSearch++;
			//DEBUG
			//for (unsigned int j=0;j<7;j++)
			//	cout << bacs[j] << " ";
			//cout << endl;	
			/*if (numBacs>3) {
				cout << "more than 3 bacs, numBacs " << numBacs << ", numPools " << numPools << endl;
				for (unsigned int j=0;j<numBacs;j++)
					cout << bacs[j] << " ";
				cout << endl;	
				//exit(1);
			};*/
			sig.setBacs(bacs, numBacs);
			pools2bacs.insert(sig);
			//DEBUG
			//if (pools2bacs.size()%1000==0)
			//	cout << "pools2bacs has " << pools2bacs.size() << " elements" << endl;
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
	for (i=0;i<numBacs;i++)
		info[k].bacs[i]=bacs[i];
	for (i=0;i<numPools;i++)
		info[k].pools[i]=pools[i];
#endif

	if (numPools<=3||(numPools>=4&&numPools<=MAXREP&&numBacs==0))
		continue;
	
	if (csets.size()==0) {//add first cset
		cset crt;
		crt.begin=k;
		crt.end=k;
		for (unsigned int i=0;i<numBacs;i++)
			crt.bacs.insert(bacs[i]);
		csets.push_back(crt);
		
#if defined DEBUG
		cout << "\t[first cset]k " << k << endl;
#endif		
		continue;
	}
	
	crtCset=csets.size() - 1;
	int crtEnd=csets[crtCset].end;

	if (numPools>MAXREP) { 
		if (k==(unsigned int)(crtEnd + 1)) { 
			csets[crtCset].end=k;
		} else {//create new cset	
			cset crt;
			crt.begin=k;
			crt.end=k;
			csets.push_back(crt);	
#if defined DEBUG
			cout << "\t[new cset]k " << k << endl;
#endif		
		}	
		
		continue;
	}	
	
	bool hasBacs=false, found=false;
	for (auto it2=csets.rbegin();it2!=csets.rend();it2++) {
		cset crt=*it2;
		
		if (crt.bacs.size()>0)
			hasBacs=true;
		
		found=false;
		for (unsigned int i=0;i<numBacs;i++) {//does this kmer belong to crt cset?
			if (crt.bacs.find(bacs[i])!=crt.bacs.end()) {
				found=true;
				break;
			}
		}

		if (found)
			break;
	}

	if (!hasBacs||found) {//if yes, merge bacs
		
		if (k==(unsigned int)(crtEnd + 1)) {//same cset	
			csets[crtCset].end=k;
		} else {//create new cset	
			cset crt;
			crt.begin=k;
			crt.end=k;
			csets.push_back(crt);	
			crtCset = crtCset + 1;
#if defined DEBUG
			cout << "\t[new cset]k " << k << endl;
#endif		
		}	
		
		for (unsigned int i=0;i<numBacs;i++)
			csets[crtCset].bacs.insert(bacs[i]);
	}

}//end for k

#if defined DEBUG
cout << endl << header << " numErrors "<< numErrors << endl;
cout<< read << endl;
for (int k=0;k<MAX_KMERS;k++) {
	cout << "kmer " << info[k].pos << " numPools " << info[k].numPools << " numBacs " << info[k].numBacs << " bacs ";
	for (i=0;i<info[k].numBacs;i++) 
		cout << info[k].bacs[i] << ",";
	//cout << endl;	
	if (info[k].numPools>=3 && info[k].numPools<=11) {
		cout << " pools ";
		for (i=0;i<info[k].numPools;i++) 
			cout << info[k].pools[i] << ",";
	}
	cout << endl;	
}
cout << header << " numErrors "<< numErrors << endl << endl;
#endif

/*use csets to fix errors by moving left and right*/

#if defined DEBUG
cout << "CORRECTION" << endl;
#endif
bool corrected=false;
bool readCorrect=true;

crtCset=0;
if ((csets.size()==(unsigned int)1)&&(csets[crtCset].begin==0)&&(csets[crtCset].end==(num-int(kmerSize))))
	crtCset=csets.size();

//read pos changed during correction
vector<int> posChanged;

//remember last k-mer "succesfully" changed for backtracking purposes
mybitsetx* lastChanged=NULL; 
//its pos in read
int lastChangedPos=-1;
//all bases we tried for last k-mer changed 
vector<char> lastChangedToTry;

while (crtCset<csets.size()) {
	//crt cset
	int bCset=csets[crtCset].begin;
	int eCset=csets[crtCset].end;
	set<unsigned short> crtCsetBacs=csets[crtCset].bacs;
	
	//next cset (or end of read) 
	//int nextECset;
	int nextBCset=num - kmerSize + 1;
	set<unsigned short> nextCsetBacs;
	if ((unsigned short)(crtCset + 1)<csets.size()) {
		nextBCset=csets[crtCset + 1].begin;
		//nextECset=csets[crtCset + 1].end;
		nextCsetBacs=csets[crtCset + 1].bacs;
	}	
	
	if (eCset==(nextBCset - 1)) {
		crtCset++;
		continue;
	}

	//prev cset (or begining of read) 
	//int prevBCset;
	int prevECset=-1;
	set<unsigned short> prevCsetBacs;
	if ((crtCset - 1)>=0) {
		//prevBCset=csets[crtCset - 1].begin;
		prevECset=csets[crtCset - 1].end;
		prevCsetBacs=csets[crtCset - 1].bacs;
	}	
	
#if defined DEBUG
	cout<< "crt cset "<< crtCset << " begin " << bCset << " end " << eCset << " bacs ";
	for (auto it=crtCsetBacs.begin();it!=crtCsetBacs.end();it++)
		cout << *it << ",";
	cout << endl; 
#endif	
	
	int direction=0;
	struct cset* neighbor=NULL;
	//if ((crtCset==0&&(bCset>0))||(crtCset>0&&((bCset - 1)!=prevECset)))
	if ((bCset - 1)!=prevECset) {
		direction=0;//go left first
		if ((crtCset - 1)>=0)
			neighbor=&csets[crtCset - 1];
	}	
	else if ((eCset + 1)!=nextBCset) {
		direction=1;//go right
		if ((unsigned short)(crtCset + 1)<csets.size())
			neighbor=&csets[crtCset + 1];

	}	
		
	int crtKmer;	
	//kmer to correct
	if (direction==0) {
		crtKmer=bCset - 1;
		e=kmers[crtKmer];
		orient=orientation[crtKmer];
	} else { 
		crtKmer=eCset + 1;
		e=kmers[crtKmer];
		orient=orientation[crtKmer];
	}	
	//kmer pools
	numPools=e.count(pools);

	mybitsetx c(e);
	
#if defined DEBUG
	cout << "\nkmer " << crtKmer << endl;
	if (posChanged.size()>0)
		cout << "Updating kmer with existing read changes" << endl;
#endif

	int pos, kmerPos;
	bool updated=false;
	//updating crtBrk kmer with possible base changes from previous fixed brk
	for (auto it2=posChanged.begin();it2!=posChanged.end();it2++) {
		pos=(*it2);
		if ((crtKmer<=pos)&&(pos<=(crtKmer+int(kmerSize)-1))) {
			updated=true;
#if defined DEBUG		
			char oldBase;
#endif			
			if (orient) {
				kmerPos=pos-crtKmer;
#if defined DEBUG		
				oldBase=c.getBase(kmerPos);
#endif			
				c.setBase(kmerPos, read[pos]);
			}	
			else {
				kmerPos=kmerSize-1-(pos-crtKmer);
#if defined DEBUG		
				oldBase=c.getBase(kmerPos);
#endif			
				c.setBase(kmerPos, mybitset::rev(read[pos]));
			}	
#if defined DEBUG
			cout<<"\t(update)crtKmer "<<crtKmer<<" read pos prev changed "<<pos<< " kmerPos updated "<<kmerPos<<" old base "<<oldBase<<" new base "<<c.getBase(kmerPos)<<endl;   
#endif
	}//end if
	}//end for
	
	//determine pos in read&kmer which needs to be changed 
	if (direction==0) {
		pos=bCset - 1;
		if (orient)
			kmerPos=0;
		else
			kmerPos=kmerSize - 1;
	} else { 
		pos=(eCset + 1) + (kmerSize - 1);
		if (orient)
			kmerPos=kmerSize - 1;
		else
			kmerPos=0;
	}	
	
	int begin, end;
	
#if defined DEBUG
	cout << "Direction " << direction << " read pos to change " << pos << " orientation " << orient << " kmerPos to change " << kmerPos << endl;
	cout << "read ";
	if (direction==0)
		begin=bCset - 1;
	else 
		begin=eCset + 1;
	for (int j=begin;j<=(begin+int(kmerSize)-1);j++)
		cout << read[j];
	cout << endl;	
	cout << "kmer " << c;
#endif	

	if (updated) {
#if defined DEBUG
		cout << "Trying " << c;
#endif		
		corrected=correct(pool, read, direction, pos, kmerPos, &csets[crtCset], neighbor, e, c, orient, pools2bacs, numSearch, numHit);
		if (corrected) {
			if (direction==0)
				{begin=bCset - 1; end=prevECset;}
			else 
				{begin=eCset + 1; end=nextBCset;}
			checkIntermediateKmers(pool, read, posChanged, direction, begin, end, &csets[crtCset], neighbor, kmers, orientation, pools2bacs, numSearch, numHit);
			continue;
		}
	}	
	
	/*try all alternatives for first/last base of this k-mer*/
	
	vector<char> toTry;
	if (crtKmer==lastChangedPos) {
		toTry=lastChangedToTry;
#if defined DEBUG
	cout << "\n!!!kmer " << crtKmer << " has been corrected before; left to try ";
	for (auto it2=toTry.begin();it2!=toTry.end();it2++) 
		cout << *it2 << " ";
	cout << endl << endl;
#endif
	} else {	
		c.flip(2*kmerPos);
		toTry.push_back(c.getBase(kmerPos));
		c.flip(2*kmerPos);
		c.flip(2*kmerPos + 1);
		toTry.push_back(c.getBase(kmerPos));
		c.flip(2*kmerPos);
		toTry.push_back(c.getBase(kmerPos));
	}

	//first alternative
	//c.flip(2*kmerPos);
	if (toTry.size()==0) {
		readCorrect=false;
		break;
	}	
	c.setBase(kmerPos, toTry.front());
	toTry.erase(toTry.begin());
	
#if defined DEBUG
	cout << "Trying " << c;
#endif	
	corrected=correct(pool, read, direction, pos, kmerPos, &csets[crtCset], neighbor, e, c, orient, pools2bacs, numSearch, numHit);
	if (corrected) {
		//record read pos which was succesfully changed
		posChanged.push_back(pos);
		
		//update last k-mer succesfully changed
		lastChanged=&c;
		lastChangedPos=crtKmer;
		lastChangedToTry=toTry;
		
		if (direction==0)
			{begin=bCset - 1; end=prevECset;}
		else 
			{begin=eCset + 1; end=nextBCset;}
		checkIntermediateKmers(pool, read, posChanged, direction, begin, end, &csets[crtCset], neighbor,  kmers, orientation, pools2bacs, numSearch, numHit);
		
		continue;
	}	
	
	//second alternative
	//c.flip(2*kmerPos);
	//c.flip(2*kmerPos+1);
	if (toTry.size()==0) {
		readCorrect=false;
		break;
	}	
	c.setBase(kmerPos, toTry.front());
	toTry.erase(toTry.begin());
	
#if defined DEBUG
	cout << "Trying " << c;
#endif	
	corrected=correct(pool, read, direction, pos, kmerPos, &csets[crtCset], neighbor, e, c, orient, pools2bacs, numSearch, numHit);
	if (corrected) {
		//record read pos which was succesfully changed
		posChanged.push_back(pos);
		
		//update last k-mer succesfully changed
		lastChanged=&c;
		lastChangedPos=crtKmer;
		lastChangedToTry=toTry;
		
		if (direction==0)
			{begin=bCset - 1; end=prevECset;}
		else 
			{begin=eCset + 1; end=nextBCset;}
		checkIntermediateKmers(pool, read, posChanged, direction, begin, end, &csets[crtCset], neighbor, kmers, orientation, pools2bacs, numSearch, numHit);
		
		continue;
	}	

	//third alternative
	//c.flip(2*kmerPos);
	if (toTry.size()==0) { 
		readCorrect=false;
		break;
	}	
	c.setBase(kmerPos, toTry.front());
	toTry.erase(toTry.begin());

#if defined DEBUG
	cout << "Trying " << c;
#endif	
	corrected=correct(pool, read, direction, pos, kmerPos, &csets[crtCset], neighbor, e, c, orient, pools2bacs, numSearch, numHit);
	if (corrected) {
		//record read pos which was succesfully changed
		posChanged.push_back(pos);
		
		//update last k-mer succesfully changed
		lastChanged=&c;
		lastChangedPos=crtKmer;
		lastChangedToTry=toTry;
		
		if (direction==0)
			{begin=bCset - 1; end=prevECset;}
		else 
			{begin=eCset + 1; end=nextBCset;}
		checkIntermediateKmers(pool, read, posChanged, direction, begin, end, &csets[crtCset],  neighbor, kmers, orientation, pools2bacs, numSearch, numHit);
	
	} else {	
		if (lastChanged) {
			if (lastChangedToTry.size()==0) {//tried all bases already
				readCorrect=false;
				break;
			}	
			
			if (direction==0)
				csets[crtCset].begin=lastChangedPos + 1;
			else 
				csets[crtCset].end=lastChangedPos - 1;
			
			posChanged.pop_back();

			continue;
		} else {
			readCorrect=false;
			break;
		}	
	}


}//end while
#if defined DEBUG
	cout << header << ":" << readCorrect << endl << endl;
#endif
out << header << ":" << readCorrect << endl << original << endl;

}//end for r
//}//end while
cout << "numSearch " << numSearch << " numHit " << numHit << endl;

out.close();
};

}

int main(int argc, char* argv[]) {

if (argc!=7)
	{
	std::cerr<<"\n\nUSE: corr <path/ReadFile> <path/MapPoolsToBacsFile> <path/OutputFile> <k-mer_size> <input_file_number> <thread>\n\n";
	exit(EXIT_FAILURE);
	}

unsigned int qgram=atoi(argv[4]);
unsigned int cores=atoi(argv[6]); 
unsigned int files=atoi(argv[5]); 

Timer* pTimer;

//load BAC sigs into trie
unsigned short* bacPools=(unsigned short*)malloc(sizeof(unsigned short)*BACS*LAYERS);
readBacMappings(bacPools);
for (unsigned short i=0;i<BACS;i++)
	trie.insert(&bacPools[i*LAYERS], 7, i+1);

class IOF ref(argv[1], argv[2], argv[3], qgram, files);
p_ref=&ref;
p_ref->ReadB();
//p_ref->CheckCorrect();

HashTable<mysig> pools2bacs;
pools2bacs.allocHash();

//int id;
struct readFasta* reads;
for (unsigned int i=36;i<39;i++) {
	
	reads=(struct readFasta*)malloc(MAX_POOL_SIZE*sizeof(struct readFasta));
	if (!reads) {
		perror("Error allocating the reads!");
		exit(1);
	}
	
	int numReads = p_ref->readReads(i, reads);
	printf("Number of reads read is %d\n", numReads);  
	
	omp_set_num_threads(cores);
	printf("Number of cores is %d\n", cores);	
	
	//correct reads in this pool
	pTimer=new Timer("correct reads in pool");
	p_ref->correctReads(i, reads, numReads, pools2bacs);
	delete pTimer;
	
	for (int r=0; r<numReads; r++) {
		free(reads[r].header);
		free(reads[r].read);
	}
	free(reads);

}//end for i

return 0;

}
