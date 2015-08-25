#include <omp.h>
#include <map>
#include <set>
#include <random>
#include <functional>

#include "trie.h"
#include "classIO.h"

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

typedef map<int, set<unsigned short>> pos2bacs;
typedef map<int, mybitsetx> pos2kmers;
typedef map<int, unsigned short> pos2orientation;

namespace IOFn {

//check intermediate k-mers between crtBrk and nextBrk -once kmer corresposnding to crtBrk was succesfully error-corrected-to make sure all k-mers are now valid and detect hidden brks 
//pos is the read pos which was succesfully changed
void IOF::checkIntermediateKmers(int pool, char* read, vector<int> posChanged, pos2bacs& breakpoints, int crtBrk, int crtBrkPools, int nextBrk, pos2kmers& kmers, pos2orientation& orientation) {

bool good=false;

unsigned short bacs[7];
unsigned short pools[POOLS];
int numBacs, numPools; 

//check all k-mers between crt brk and the next one
//because there might be clustered errors which obfuscated some brks
for (int p=crtBrk+1;p<nextBrk;p++) {
	mybitsetx e=kmers[p];
	unsigned short orient=orientation[p];
	
	//update pool counts in hashtable
	sh.searchUpdatePool(e, pool, -1);
	
	//DEBUG
	//cout << "[checkIntermediateKmers]Updating kmer with existing read changes" << endl;
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
			//DEBUG
			//cout<<"\t[checkIntermediateKmers]kmer p "<<p<<" read pos prev changed "<<pos<<" new read base "<<read[pos]<< " kmerPos changed "<<kmerPos<<" old kmer base "<<oldBase<<" new kmer base "<<e.getBase(kmerPos)<<" orientation "<<orient<<endl;   
		}//end if
	}//end for it	
	
	mybitsetx g;
	mybitsetx o(e);
	e.invert();
	if (o<e)
		g=o;
	else
		g=e;

	//update pool counts in hashtable
	sh.searchUpdatePool(g, pool, 1);
	
	good=false;
	numBacs=0; numPools=0;
	if (sh.searchCopy(g)) {//search kmer in hashtable
		numPools=g.count(pools);
		
		if (numPools>25/*&&crtBrkPools>25*/) 
			good=true;
		
		if (numPools>=4&&numPools<=25) {
			trie.search(pools, numPools, bacs,numBacs);//search bacs in trie
			for (auto it=breakpoints.begin();it!=breakpoints.end();it++) {
				int brk=it->first;
				set<unsigned short> brkBacs=it->second;
				for (int i=0;i<numBacs;i++) {
					if (brkBacs.find(bacs[i])!=brkBacs.end()) {
						good=true;
						break;
					}
				}
				
				if (good)
					break;
			}
		}	
	}
	
	//if (good) {
		//DEBUG
		//cout << "\t[checkIntermediateKmers]checked kmer at pos  " << p << endl; 
		//cout << "\t[checkIntermediateKmers]numPools " << numPools << " crtBrkPools " << crtBrkPools << endl;
		//cout << "\t[checkIntermediateKmers]numBacs " << numBacs << ": ";
		//for (int i=0;i<numBacs;i++) 
		//	cout << bacs[i] << ",";
		//cout << endl;	
	//}

	if (!good) {//a potential hidden brk
		//DEBUG
		cout << "\t[checkIntermediateKmers]checked kmer at pos  " << p << endl; 
		cout << "\t[checkIntermediateKmers]numPools " << numPools << " crtBrkPools " << crtBrkPools << endl;
		cout << "\t[checkIntermediateKmers]numBacs " << numBacs << ": ";
		for (int i=0;i<numBacs;i++) 
			cout << bacs[i] << ",";
		cout << endl;	
		cout << "\t[checkIntermediateKmers]new breakpoint added at " << p << endl;
		breakpoints.insert(pair<int,set<unsigned short>>(p,set<unsigned short>(bacs, bacs+numBacs)));
		break;
	}
}//end for p

};


bool IOF::correct(int pool, char* read, int pos, pos2bacs& breakpoints, int crtBrk, int& crtBrkPools, const mybitsetx& e, mybitsetx c, unsigned short orient) {

int i;
bool correct=false;
unsigned short bacs[7];
unsigned short pools[POOLS];
int numBacs, numPools; 

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
	//DEBUG
	cout << "kmer not found in the hashtable!" << endl;
	return false;
}	

crtBrkPools=numPools=g.count(pools);
if (numPools<4) {
	//DEBUG
	cout << "[correct]no pools found " << numPools << endl;
	return false;	
}

if (numPools>25)
	correct=true;

numBacs=0;
if (numPools>=4&&numPools<=25) {
	trie.search(pools, numPools, bacs, numBacs);
	if (numBacs==0)	{
		//DEBUG
		cout << "[correct]no bacs found " << numBacs << endl;
		return false;
	}
	
	set<unsigned short> crtBrkBacs;
	for (i=0;i<numBacs;i++) 
		crtBrkBacs.insert(bacs[i]);
	
	int brk;
	set<unsigned short> brkBacs;
	//check that bacs found agree with at least another breakpoint
	for (auto it2=breakpoints.begin();it2!=breakpoints.end();it2++) {
		brk =it2->first;
		if (brk==crtBrk)
			continue;
		brkBacs=it2->second;
		
		bool found=false;
		for (i=0;i<numBacs;i++) {
			if (brkBacs.find(bacs[i])!=brkBacs.end()) {	
				found=true;
				break;
			}
		}
			
		if (found) {
			//merge the two breakpoints
			crtBrkBacs.insert(brkBacs.begin(), brkBacs.end());
			//delete one breakpoint
			breakpoints.erase(brk);
			correct=true;
		}	
	}//end for	
	
	breakpoints[crtBrk]=crtBrkBacs;

}//if (numPools>=4&&numPools<=25) {

//DEBUG
cout << "[correct]Found " << g; 
cout << "[correct]numPools " << numPools << endl;
cout << "[correct]numBacs " << numBacs << ": ";
for (i=0;i<numBacs;i++) 
	cout << bacs[i] << ",";
cout << endl;	

if (correct) {
	//change base in read
	char newBase; 
	char oldBase=read[pos];
	int kmerPos=pos-crtBrk;
	if (orient) 
		newBase=o.getBase(kmerPos);
	else 
		newBase=c.getBase(kmerPos);
	read[pos]=newBase;
	//DEBUG
	cout << "[correct]read pos " << pos << " old base " << oldBase << " new base " << newBase << endl;
	//update pool counts in hashtable
	sh.searchUpdatePool(e, pool, -1);
	sh.searchUpdatePool(g, pool, 1);
	
	return true;
}
	
return false;

};


void IOF::correctReads(int pool, struct readFasta* reads, int numReads)
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

int i;
struct kmerInfo info[MAX_KMERS];

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
	
	pos2kmers kmers;
	pos2bacs breakpoints;
	pos2orientation orientation;
	unsigned short bacs[7];
	unsigned short pools[POOLS];
	unsigned short orient;
	int numBacs, numPools; 
	unsigned short haveBacs=0;
	int crtBrk=-1, crtBrkPools=0;
	set<unsigned short> crtBrkBacs;

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
	
	kmers.insert(pair<int, mybitsetx>(k, e));
	orientation.insert(pair<int, unsigned short>(k, orient));
	
	numBacs=0;
	//get positive pools and their number for this k-mer
	numPools=e.count(pools);

	if (numPools>=4 && numPools<=25)//get bacs for this kmer 
		trie.search(pools, numPools, bacs, numBacs);

	//DEBUG
	info[k].pos=k;
	info[k].numPools=numPools;
	info[k].numBacs=numBacs;
	for (i=0;i<numBacs;i++)
		info[k].bacs[i]=bacs[i];
	for (i=0;i<numPools;i++)
		info[k].pools[i]=pools[i];
			
	if (crtBrk!=-1) 
		crtBrkBacs = breakpoints[crtBrk];
	
	if (numPools<=3) {//this unique k-mer is most probably an error
		if ((crtBrk==-1) || (crtBrkPools>3) ||(crtBrkBacs.size()>0)) {
			crtBrk=k;//create new breakpoint
			crtBrkPools=numPools;
			breakpoints.insert(pair<int, set<unsigned short>>(crtBrk, set<unsigned short>()));
						
			//DEBUG
			//if (numErrors==0&&k>0) 
				cout << "\t[sig<=3]crtBrk " << crtBrk << endl;
		}
		continue;
	}	

	if (crtBrkBacs.empty()) {//the crt brk point is associated with a unique k-mer
		if ((crtBrk==-1) || (numBacs>0)) { 
			crtBrk=k;//create new breakpoint
			crtBrkPools=numPools;
			breakpoints.insert(pair<int, set<unsigned short>>(crtBrk, set<unsigned short>(bacs, bacs + numBacs)));
			if (numBacs>0)
				haveBacs++;
			//DEBUG
			//if (numErrors==0&&k>0) 
				cout << "\t[empty crtBrkBacs]crtBrk " << crtBrk << endl;
		}
		continue;
	}	
	
	if (numPools>25)
		continue;
	
	/*if (numBacs==0) {
		crtBrk=k;
		crtBrkPools=numPools;
		kmers.insert(pair<int, mybitsetx>(crtBrk, e));
		orientation.insert(pair<int, unsigned short>(crtBrk, orient));
		breakpoints.insert(pair<int, set<unsigned short>>(crtBrk, set<unsigned short>()));
		//DEBUG
		//if (numErrors==0&&k>0) 
			cout << "\t[numBacs=0]crtBrk " << crtBrk << endl;
		continue;
	}*/

	bool found=false;//crt brk point not empty
	//check if crt Bacs can be merged to any of the prev breakpoint bacs
	//for (auto it=breakpoints.begin();it!=breakpoints.end();it++) {
	//	int brk=it->first;
	//	set<unsigned short> brkBacs=it->second;
	//	if (brkBacs.size()==0)
	//		continue;
	for (int i=0;i<numBacs;i++) {
		if (crtBrkBacs.find(bacs[i])!=crtBrkBacs.end()) {
			found=true;
			break;
		}
	}	
	
	if (found) {//merge
		for (int i=0;i<numBacs;i++)
			breakpoints[crtBrk].insert(bacs[i]);
	}
	//}
	
	if (!found) {//create new breakpoint
		crtBrk=k;
		crtBrkPools=numPools;
		if (numBacs>0) {
			breakpoints.insert(pair<int, set<unsigned short>>(crtBrk, set<unsigned short>(bacs, bacs + numBacs)));
			haveBacs++;
		}
		else {
			breakpoints.insert(pair<int, set<unsigned short>>(crtBrk, set<unsigned short>()));
		}
		//DEBUG
		//if (numErrors==0&&k>0) 
			cout << "\t[conflicting crtBrkBacs]crtBrk " << crtBrk << endl;
	}

}//end for k

int numBrk=breakpoints.size();
if ((numBrk==1) && (breakpoints.find(0)!=breakpoints.end())) 
	numBrk=0;

//DEBUG
//if (numErrors==0&&numBrk>0) {

cout << endl << header << " numErrors "<< numErrors << " numBrk " << numBrk << endl;
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
cout << header << " numErrors "<< numErrors << " numBrk " << numBrk << endl << endl;

//}//end if

/*iterate over and fix breakpoints*/

cout << "CORRECTION" << endl;
bool corrected=false;
bool readCorrect=true;
vector<int> posChanged;
for (auto it=breakpoints.begin();it!=breakpoints.end();it++) {
	
	//crt brk and associated bacs
	int crtBrk =it->first;
	crtBrkBacs=it->second;
	//kmer associated to crtBrk
	e=kmers[crtBrk];
	//kmer orientation
	orient=orientation[crtBrk];
	//kmer pool
	numPools=e.count(pools);
	
	if ((crtBrk==0)&&(numBrk==0||crtBrkBacs.size()>0||numPools>25))
		continue;
	if ((haveBacs==1)&&(crtBrkBacs.size()>0))
		continue;

	//next brk (or end of read) and associated bacs
	int nextBrk=num-kmerSize+1;
	set<unsigned short> nextBrkBacs;
	it++;
	if (it!=breakpoints.end()) {
		nextBrk=it->first;
		nextBrkBacs=it->second;
	}	
	it--;
	
	//DEBUG
	cout << "crtBrk " << crtBrk << " haveBacs " << haveBacs << " bacs ";
	for (auto it=crtBrkBacs.begin();it!=crtBrkBacs.end();it++)
		cout << *it << ",";
	cout << endl; 
	cout << "read "; 
	for (int j=crtBrk;j<=(crtBrk+int(kmerSize)-1);j++)
		cout << read[j];
	cout << endl;	
	cout << "kmer " << e;
	
	mybitsetx c(e);
	
	//DEBUG
	if (posChanged.size()>0)
		cout << "Updating kmer with existing read changes" << endl;
	
	int pos, kmerPos;
	bool updated=false;
	//updating crtBrk kmer with possible base changes from previous fixed brk
	for (auto it2=posChanged.begin();it2!=posChanged.end();it2++) {
		pos=(*it2);
		if ((crtBrk<=pos)&&(pos<=(crtBrk+int(kmerSize)-1))) {
			updated=true;
			char oldBase;
			if (orient) {
				kmerPos=pos-crtBrk;
				oldBase=c.getBase(kmerPos);
				c.setBase(kmerPos, read[pos]);
			}	
			else {
				kmerPos=kmerSize-1-(pos-crtBrk);
				oldBase=c.getBase(kmerPos);
				c.setBase(kmerPos, mybitset::rev(read[pos]));
			}	
			//DEBUG
			cout<<"\tcrtBrk "<<crtBrk<<" read pos prev changed "<<pos<< " kmerPos changed "<<kmerPos<<" old base "<<oldBase<<" new base "<<c.getBase(kmerPos)<<endl;   
		}//end if
	}//end for	
	
	
	/*try all alternatives for first/last base of this k-mer*/
	
	//determine pos in read which needs to be changed 
	pos=crtBrk + int(kmerSize) -1;
	if (crtBrk<int(kmerSize)) {//deal with special case when err in first kmerSize bases of read
		if ((nextBrk<int(kmerSize))&&(nextBrkBacs.size()>0))
			pos=nextBrk-1;
	}
	//determine pos in k-mer which needs to be changed
	if (orient)
		kmerPos=pos-crtBrk;
	else
		kmerPos=kmerSize-1-(pos-crtBrk);
	//DEBUG
	cout << "\tpos to change " << pos << " kmerPos to change " << kmerPos << endl;

	if (updated) {
		//DEBUG
		cout << "Trying " << c;
		corrected=correct(pool, read, pos, breakpoints,crtBrk,numPools,e,c,orient);
		if (corrected) {
			checkIntermediateKmers(pool, read, posChanged, breakpoints, crtBrk, numPools, nextBrk, kmers, orientation);
			continue;
		}
	}	
	
	//first alternative
	c.flip(2*kmerPos);
	//DEBUG
	cout << "Trying " << c;
	corrected=correct(pool, read, pos, breakpoints,crtBrk,numPools,e,c,orient);
	if (corrected) {
		//record read pos which was succesfully changed
		posChanged.push_back(pos);
		//this may discover hidden brks and therefore modify the breakpoints map  
		checkIntermediateKmers(pool, read, posChanged, breakpoints, crtBrk, numPools, nextBrk, kmers, orientation);
		continue;
	}	
	
	//second alternative
	c.flip(2*kmerPos);
	c.flip(2*kmerPos+1);
	//DEBUG
	cout << "Trying " << c;
	corrected=correct(pool, read, pos, breakpoints,crtBrk,numPools,e,c,orient);
	if (corrected) {
		//record read pos which was succesfully changed
		posChanged.push_back(pos);
		//this may discover hidden brks and therefore modify the breakpoints map  
		checkIntermediateKmers(pool, read, posChanged, breakpoints, crtBrk, numPools, nextBrk, kmers, orientation);
		continue;
	}	

	//third alternative
	c.flip(2*kmerPos);
	//DEBUG
	cout << "Trying " << c;
	corrected=correct(pool, read, pos, breakpoints,crtBrk,numPools,e,c,orient);
	if (corrected) {
		//record read pos which was succesfully changed
		posChanged.push_back(pos);
		//this may discover hidden brks and therefore modify the breakpoints map  
		checkIntermediateKmers(pool, read, posChanged, breakpoints, crtBrk, numPools, nextBrk, kmers, orientation);
	} else {	
		readCorrect=false;
		break;
	}


}//end for
//DEBUG
cout << header << ":" << readCorrect << endl << endl;
out << header << ":" << readCorrect << endl << original << endl;

}//end for r
//}//end while

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

//load BAC sigs into trie
unsigned short* bacPools=(unsigned short*)malloc(sizeof(unsigned short)*BACS*LAYERS);
readBacMappings(bacPools);
for (unsigned short i=0;i<BACS;i++)
	trie.insert(&bacPools[i*LAYERS], 7, i+1);

class IOF ref(argv[1], argv[2], argv[3], qgram, files);
p_ref=&ref;
p_ref->ReadB();

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
	p_ref->correctReads(i, reads, numReads);
	
	for (int r=0; r<numReads; r++) {
		free(reads[r].header);
		free(reads[r].read);
	}
	free(reads);

}//end for i

return 0;

}
