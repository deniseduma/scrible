#include <random>
#include <functional>

#include "Timer.h"
#include "trie.h"
#include "classIO.h"

//#define DEBUG 1

//Best so far.

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
void IOF::checkIntermediateKmers(int pool, char* read, const vector<int>& posChanged, int direction, int begin, int end, struct cset& crt, const struct cset& neighbor, const vector<mybitsetx>& kmers, const vector<unsigned short>& orientation, HashTable<mysig>& pools2bacs, unsigned int& numSearch, unsigned int& numHit) {}

bool IOF::correct(int pool, char* read, int direction, int pos, int kmerPos, unsigned short crtCset, vector<struct cset>& csets, const mybitsetx& e, mybitsetx c, unsigned short orient, HashTable<mysig>& pools2bacs, unsigned int& numSearch, unsigned int& numHit, unsigned short* b, unsigned int& nB) {

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
		cout << "[correct]numPools " << numPools << ": ";
		for (unsigned int i=0;i<numPools;i++)
			cout << pools[i] << " ";
		cout << endl;	
#endif		
		return false;
	}
	
	//check that bacs found agree with other cset bacs
	bool hasBacs=false, found=false;
	for (auto it=csets.begin();it!=csets.end();it++) {
		if (it->bacs.size()==0)
			continue;
		
		hasBacs=true;
		
		found=false;
		for (unsigned int i=0;i<numBacs;i++) {//does this kmer belong to crt cset?
			if (it->bacs.find(bacs[i])!=it->bacs.end()) {
				found=true;
				break;
			}
		}

		if (found)
			break;
	};

	if (!hasBacs||found) {
		correct=true;
		for (i=0;i<numBacs;i++) 
			csets[crtCset].bacs.insert(bacs[i]);
	}	
	
}//if (numPools>=LOW&&numPools<=HIGH) {

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
		csets[crtCset].begin--;
	else	
		csets[crtCset].end++;
	//change base in read
	char newBase;
	char oldBase=read[pos];
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
	cout << "[correct]crtCset begin " << csets[crtCset].begin << " crtCset end " << csets[crtCset].end << endl; 
#endif
	return true;
}//if (correct)
	
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

//remember last k-mer "succesfully" changed for backtracking purposes
mybitsetx* lastChanged; 
//its pos in read
int lastChangedPos;

for (int r=0;r<numReads;r++) { 
	
	char* header = reads[r].header;
	char* read = reads[r].read;
	unsigned int num = strlen(read);
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
	unsigned int numBacs;
	unsigned short bacs[MAX_BACS];
	unsigned short orient;
	vector<unsigned short> orientation;

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
	
	bool inHash=sh.searchCopy(e); 
	
	kmers.push_back(e);
	orientation.push_back(orient);
	
	//get positive pools and their number for this k-mer
	numPools=e.count(pools);

	//unsigned short bacs[3];
	//for (unsigned int j=0;j<7;j++)
	//	cout << bacs[j] << " ";
	//cout << endl;	
	numBacs=0;
	if (numPools>=LOW && numPools<=HIGH) { //get bacs for this kmer 
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

	if (numPools<=3||(numPools>=LOW&&numPools<=HIGH&&numBacs==0))
		continue;
	
	if (csets.size()==0) {//add first cset
		cset crt;
		crt.begin=k;
		crt.end=k;
		for (unsigned int i=0;i<numBacs;i++)
			crt.bacs.insert(bacs[i]);
		csets.push_back(crt);
		
#if defined DEBUG
		cout << "\t[first cset]at k " << k << endl;
#endif		
		continue;
	}
	
	crtCset=csets.size() - 1;
	int crtEnd=csets[crtCset].end;

	if (numPools>HIGH) { 
		if (k==(unsigned int)(crtEnd + 1)&&csets[crtCset].bacs.size()==0) { 
			csets[crtCset].end=k;
		} else {//create new cset	
			cset crt;
			crt.begin=k;
			crt.end=k;
			csets.push_back(crt);	
#if defined DEBUG
			cout << "\t[new cset] at k " << k << endl;
#endif		
		}	
		
		continue;
	}	
	
	bool hasBacs=false, found=false;
	for (auto it2=csets.rbegin();it2!=csets.rend();it2++) {
		cset crt=*it2;
		
		if (crt.bacs.size()==0)
			continue;
		
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
		
		if (k==(unsigned int)(crtEnd + 1)&&csets[crtCset].bacs.size()>0) {//same cset	
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
	/*if (info[k].numPools>=3 && info[k].numPools<=11) {
		cout << " pools ";
		for (i=0;i<info[k].numPools;i++) 
			cout << info[k].pools[i] << ",";
	}*/
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
if ((csets.size()==(unsigned int)1)&&(csets[crtCset].begin==0)&&((unsigned int)csets[crtCset].end==(num-kmerSize)))
	crtCset=csets.size();

int crtKmer;
//read pos changed during correction
vector<int> posChanged;

while (crtCset<csets.size()) {
	//crt cset
	int bCset=csets[crtCset].begin;
	int eCset=csets[crtCset].end;
	
	//prev cset or begining of read 
	//int prevBCset;
	int prevECset=-1;
	//set<unsigned short> prevCsetBacs;
	if ((crtCset - 1)>=0) {
		//prevBCset=csets[crtCset - 1].begin;
		prevECset=csets[crtCset - 1].end;
		//prevCsetBacs=csets[crtCset - 1].bacs;
	}	
	
	//next cset or end of read 
	int nextBCset=num - kmerSize + 1;
	//int nextECset;
	//set<unsigned short> nextCsetBacs;
	if ((unsigned short)(crtCset + 1)<csets.size()) {
		nextBCset=csets[crtCset + 1].begin;
		//nextECset=csets[crtCset + 1].end;
		//nextCsetBacs=csets[crtCset + 1].bacs;
	}	
	
	if (((bCset - 1)==prevECset)&&((eCset + 1)==nextBCset)) {
		crtCset++;
		continue;
	}


#if defined DEBUG
	set<unsigned short> crtCsetBacs=csets[crtCset].bacs;
	cout<< "[CORRECTION]crt cset "<< crtCset << " begin " << bCset << " end " << eCset << " bacs ";
	for (auto it=crtCsetBacs.begin();it!=crtCsetBacs.end();it++)
		cout << *it << ",";
	cout << endl; 
#endif	
	
	//Compute direction left or right
	int direction=0;
	if ((bCset - 1)!=prevECset) 
		direction=0;//go left first
	else if ((eCset + 1)!=nextBCset) 
		direction=1;//go right
	
	lastChanged=NULL; 
	lastChangedPos=-1;
	vector<char> lastChangedToTry;
	while (true) {
		
		if (direction==0) {
			assert(csets[crtCset].begin>prevECset);
			crtKmer=csets[crtCset].begin - 1;
			if (crtKmer==prevECset)
				//direction=1;
				break;
			e=kmers[crtKmer];
			orient=orientation[crtKmer];
		}	
		if (direction==1) { 
			assert(csets[crtCset].end<nextBCset);
			crtKmer=csets[crtCset].end + 1;
			if (crtKmer==nextBCset) 
				break;
			e=kmers[crtKmer];
			orient=orientation[crtKmer];
		}	
		//kmer pools
		//numPools=e.count(pools);

		mybitsetx c(e);//keep a copy of the original k-mer
	
#if defined DEBUG
	cout << "\nkmer " << crtKmer << endl;
	if (posChanged.size()>0)
		cout << "Updating kmer with existing read changes" << endl;
#endif

		int pos, kmerPos;
		//bool updated=false;
		//updating crt kmer with possible base changes from previous fixed kmers
		for (auto it2=posChanged.begin();it2!=posChanged.end();it2++) {
			assert(pos>=0&&pos<(int)num);
			/*if (!(pos>=0&&pos<(int)num)) {
				cout << header << endl;
				exit(1);
			}*/
			pos=(*it2);
			if ((crtKmer<=pos)&&(pos<=(crtKmer+int(kmerSize)-1))) {
				//updated=true;
#if defined DEBUG		
				char oldBase;
#endif			
				if (orient) {
					kmerPos=pos-crtKmer;
#if defined DEBUG		
					oldBase=c.getBase(kmerPos);
#endif			
					c.setBase(kmerPos, read[pos]);
				} else {
					kmerPos=kmerSize-1-(pos-crtKmer);
#if defined DEBUG		
					oldBase=c.getBase(kmerPos);
#endif			
					c.setBase(kmerPos, mybitset::rev(read[pos]));
				}	
#if defined DEBUG
	cout<<"\t(update)crtKmer "<<crtKmer<<" read pos prev changed "<<pos<< " kmerPos updated "<<kmerPos<<" old base "<<oldBase<<" new base "<<c.getBase(kmerPos)<<endl<<endl;   
#endif
			}//end if
		}//end for
	
		//determine pos in read&kmer which needs to be changed 
		if (direction==0) {
			pos=csets[crtCset].begin - 1;
			assert(pos>prevECset);
			if (orient)
				kmerPos=0;
			else
				kmerPos=kmerSize - 1;
		} else { 
			pos=csets[crtCset].end + kmerSize;
			assert((csets[crtCset].end + 1)<nextBCset);
			if (orient)
				kmerPos=kmerSize - 1;
			else
				kmerPos=0;
		}	
	
#if defined DEBUG
	cout << "Direction " << direction << " read pos to change " << pos << " orientation " << orient << " kmerPos to change " << kmerPos << endl;
	int begin;
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

#if defined DEBUG
	cout << "Trying same k-mer " << c;
#endif		
		//if (updated) {
		corrected=correct(pool, read, direction, pos, kmerPos, crtCset, csets, e, c, orient, pools2bacs, numSearch, numHit, bacs, numBacs);
		if (corrected) 
			continue;
		//}	
		
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
		corrected=correct(pool, read, direction, pos, kmerPos, crtCset, csets, e, c, orient, pools2bacs, numSearch, numHit, bacs, numBacs);
		if (corrected) {
			//record read pos which was succesfully changed
			posChanged.push_back(pos);
			//assert(posChanged.size()<=5);
			if (posChanged.size()>5) {
				readCorrect=false;
				break;
			}
		
			//update last k-mer succesfully changed
			lastChanged=&c;
			lastChangedPos=crtKmer;
			lastChangedToTry=toTry;
			continue;
		}	
	
		if (toTry.size()==0) {
			readCorrect=false;
			break;
		}	
		c.setBase(kmerPos, toTry.front());
		toTry.erase(toTry.begin());
	
#if defined DEBUG
	cout << "Trying " << c;
#endif	
		corrected=correct(pool, read, direction, pos, kmerPos, crtCset, csets, e, c, orient, pools2bacs, numSearch, numHit, bacs, numBacs);
		if (corrected) {
			//record read pos which was succesfully changed
			posChanged.push_back(pos);
			//assert(posChanged.size()<=5);
			if (posChanged.size()>5) {
				readCorrect=false;
				break;
			}
		
			//update last k-mer succesfully changed
			lastChanged=&c;
			lastChangedPos=crtKmer;
			lastChangedToTry=toTry;
			continue;
		}	

		if (toTry.size()==0) { 
			readCorrect=false;
			break;
		}	
		c.setBase(kmerPos, toTry.front());
		toTry.erase(toTry.begin());

#if defined DEBUG
	cout << "Trying " << c;
#endif	
		corrected=correct(pool, read, direction, pos, kmerPos, crtCset, csets, e, c, orient, pools2bacs, numSearch, numHit, bacs, numBacs);
		if (corrected) {
			//record read pos which was succesfully changed
			posChanged.push_back(pos);
			//assert(posChanged.size()<=5);
			if (posChanged.size()>5) {
				readCorrect=false;
				break;
			}
		
			//update last k-mer succesfully changed
			lastChanged=&c;
			lastChangedPos=crtKmer;
			lastChangedToTry=toTry;
		
	
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
				
				//revert pool counts in hashtable
				c=*lastChanged;
				mybitsetx o(c);
				c.invert();
				mybitsetx g;
				if (c<o)
					g=c;
				else
					g=o;
				sh.searchUpdatePool(g, pool, -1);
				sh.searchUpdatePool(kmers[lastChangedPos], pool, 1);
				
				continue;
			} else {
				readCorrect=false;
				break;
			}	
		}//end else

	}//end while(true)

	if (readCorrect==false)
		break;
	
	assert((csets[crtCset].begin - 1)==prevECset);	
	//assert((csets[crtCset].end + 1)==nextBCset);
	
	//crtCset++;
	
	/*if ((csets[crtCset].begin - 1)!=prevECset) {
		cout << "bCset " << csets[crtCset].begin << ", prevECset " << prevECset << endl;
		cout << header << endl;
		//exit(1);
	}
	if ((csets[crtCset].end + 1)!=nextBCset) {
		cout << "eCset " << csets[crtCset].end << ", nextBCset " << nextBCset << endl;
		cout << header << endl;
		//exit(1);
	}*/

}//end while
#if defined DEBUG
	cout << header << ":" << readCorrect << endl << endl;
#endif
out << header << ":" << readCorrect << endl << read << endl;

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
for (unsigned int i=36;i<37;i++) {
	
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
