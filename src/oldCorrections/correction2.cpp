#include <omp.h>
#include <map>
#include <set>
#include <random>
#include <functional>

#include "trie.h"
#include "classIO.h"

#define MAX_POOL_SIZE 2000000

unsigned char c[8]={ 0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80};
unsigned char r[8]={ 0xFE, 0xFD, 0xFB, 0xF7, 0xEF, 0xDF, 0xBF, 0x7F};

unsigned int hashSizes[SIZES]={1000393, 2500009, 5000549, 10000643, 20000003, 30000001, 40000003, 50000017, 60000011,70000027, 80000023, 90000049, 100000567, 150000001, 200000477, 250000013, 300000641, 350000041, 400000651, 450000007, 500000461, 550000001, 600000607, 650000011, 700000747, 750000007, 800000533, 850000013, 900000637, 950000017, 1000000531, 1250000027, 1500000713, 1750000027, 2000000579, 2147483647}; 

using namespace std;
using namespace general;
using namespace IOFn;

Trie trie;
class IOF* p_ref;

typedef map<int, set<unsigned short>> mymap;

namespace IOFn {

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
	char* errs=strchr(header, ':');
	int numErrors=0;
	//count the number of errors the read contains
	while(*errs) if (*errs++ == '_') numErrors++;
	int num = strlen(read);
	
	mymap breakpoints;
	int numBacs=0, crtBrk=-1, crtBrkPools=0;
	unsigned short bacs[7];
	unsigned short pools[POOLS];
	set<unsigned short> crtBrkBacs;

	class mybitsetx b,d,e,g;
	for (unsigned int k=0; k<=num-kmerSize+1; k++) {
		int z=0;
		int j=kmerSize*2-1;
		bool findN=false;
		for (unsigned int c=0; c<kmerSize; c++) {
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
			default:
			findN=true;
			c=kmerSize;
			}
	}//end for c		
		
	if (findN) 
		continue;
	if (d<b) 
		e=d;
	else
		e=b;
	
	bool inHash=sh.searchCopy(e); 
	
	numBacs=0;
	//get positive pools and their number for this k-mer
	int numPools=e.count(pools);

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
			if (numErrors==0&&k>0) 
				cout << "\t[sig<=3]crtBrk " << crtBrk << endl;
		}
		continue;
	}	

	if (crtBrkBacs.empty()) {//the crt brk point is associated with a unique k-mer
		if ((crtBrk==-1) || (numBacs>0)) { 
			crtBrk=k;//create new breakpoint
			crtBrkPools=numPools;
			breakpoints.insert(pair<int, set<unsigned short>>(crtBrk, set<unsigned short>(bacs, bacs + numBacs)));
			//DEBUG
			if (numErrors==0&&k>0) 
				cout << "\t[empty crtBrkBacs]crtBrk " << crtBrk << endl;
		}
		continue;
	}	
	
	if (numPools>25)
		continue;

	bool found=false;//crt brk point not empty
	for (int i=0;i<numBacs;i++) {
		if (crtBrkBacs.find(bacs[i])!=crtBrkBacs.end()) {
			found=true;
			break;
		}
	}	
	if (found) {
		for (int i=0;i<numBacs;i++)
			breakpoints[crtBrk].insert(bacs[i]);
	} else {
		crtBrk=k;
		crtBrkPools=numPools;
		if (numBacs>0)
			breakpoints.insert(pair<int, set<unsigned short>>(crtBrk, set<unsigned short>(bacs, bacs + numBacs)));
		else
			breakpoints.insert(pair<int, set<unsigned short>>(crtBrk, set<unsigned short>()));
			//DEBUG
			if (numErrors==0&&k>0) 
				cout << "\t[conflicting crtBrkBacs]crtBrk " << crtBrk << endl;
	}

}//end for k
//int numBrk=(breakpoints.size()>0) ? (breakpoints.size() - 1) : 0;
int numBrk=breakpoints.size();
if ((numBrk==1) && (breakpoints.find(0)!=breakpoints.end())) 
	numBrk=0;
out<<header<<":"<<numBrk<<endl<<read<<endl;

//DEBUG
if (numErrors==0&&numBrk>0) {

cout << endl << header << " numErrors "<< numErrors << endl;
for (int k=0;k<MAX_KMERS;k++) {
	cout << "kmer " << info[k].pos << " numPools " << info[k].numPools << " numBacs " << info[k].numBacs << " bacs ";
	for (i=0;i<info[k].numBacs;i++) 
		cout << info[k].bacs[i] << ",";
	cout << endl;	
	if (info[k].numPools>=3 && info[k].numPools<=12) {
		cout << "(pools) ";
		for (i=0;i<info[k].numPools;i++) 
			cout << info[k].pools[i] << ",";
		cout << endl;	
	}
}
cout<<header<<":"<<numBrk<<endl<<read<<endl<<endl;

}//end if

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

int id, r;
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
	
	for (r=0; r<numReads; r++) {
		free(reads[r].header);
		free(reads[r].read);
	}
	free(reads);

}//end for i

return 0;

}
