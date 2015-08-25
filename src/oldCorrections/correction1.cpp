#include <omp.h>
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

namespace IOFn {

void IOF::correctReads(int pool, struct readFasta* reads, int numReads)
{

//output files for corrected and uncorrected reads
ostringstream of, of1;

char delimC[] = ".";
Parser parser;
parser.update(delimC, ReadFname);

ofstream out, out1;
of<<parser.get(0)<<pool<<"."<<"corr";
out.open(of.str().c_str(), ios_base::out);
if(!out) {
	cerr<<"Error opening output file "<<of.str()<<endl;
	exit(EXIT_FAILURE);
}	
of1<<parser.get(0)<<pool<<"."<<"err";
out1.open(of.str().c_str(), ios_base::out);
if(!out1) {
	cerr<<"Error opening output file "<<of1.str()<<endl;
	exit(EXIT_FAILURE);
}	

bool readChanged=true;
//while (readChanged) {
readChanged=false;
for (int r=0;r<numReads;r++) { 
	
	char* header = reads[r].header;
	char* read = reads[r].read;
	int num = strlen(read);
	//DEBUG
	printf("%s\n", header);
	
	//int numKmers=0;
	//int kmers[POOLS*MAX_KMERS];
	
	unsigned short bacs[7];
	unsigned short cBacs[7]; 
	unsigned short numBacs, numCBacs, maxDist;

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
	if (!sh.searchCopy(e)) 
		continue;
	
	numBacs=0;
	int numPools=e.count();
	if (numPools()>=6) //search bacs for this kmer 
		trie.search(e.getPools(), numPools, bacs, numBacs);
		
	bool hasBacs = (numBacs>0) ? true : false;
	
	if (hasBacs && numCBacs==0) {
		for (i=0;i<numBacs;i++)
			cBacs[numCBacs++]=bacs[i];
		continue;
	}		
	
	bool found=false;
	for (i=0;i<numBacs;i++) {
		for (j=0;j<numCBacs;j++) {
			if (bacs[i]==cBacs[j]) {
				found=true;
				break;
			}
		}	
	}
	
	if (!hasBacs || !found) {
		//try all 3 posibilities for first base of kmer
		mybitsetx c=e; 
			
		c.flip(0); 
		if (sh.searchCopy(c)) {
			trie.search(c.getPools(), c.count(), 5, bacs, numBacs);
			for (i=0;i<numBacs;i++) 
				for (j=0;j<numCBacs;j++)
					if (bacs[i]==cBacs[j]) { 
						//change base in read
						read[k]=c.getBase(0);
						readChanged=true;
						//update pool counts in hashtable
						sh.searchUpdatePool(e, pool, -1);
						sh.searchUpdatePool(c, pool, 1);
						//continue to next k-mer
						continue;
					}		
		}	
			
		c.flip(0);
		c.flip(1);
		if (sh.searchCopy(c)) {
			trie.search(c.getPools(), c.count(), 5, bacs, numBacs);
			for (i=0;i<numBacs;i++) 
				if (bacs[i]==bacSoFar) { 
					//change base in read
					read[k]=c.getBase(0);
					readChanged=true;
					//update pool counts in hashtable
					sh.searchUpdatePool(e, pool, -1);
					sh.searchUpdatePool(c, pool, 1);
					//continue to next k-mer
					continue;
				}		
		}	
			
		c.flip(0);
		if (sh.searchCopy(c)) {
			trie.search(c.getPools(), c.count(), 5, bacs, numBacs);
			for (i=0;i<numBacs;i++) 
				if (bacs[i]==bacSoFar) { 
					//change base in read
					read[k]=c.getBase(0);
					readChanged=true;
					//update pool counts in hashtable
					sh.searchUpdatePool(e, pool, -1);
					sh.searchUpdatePool(c, pool, 1);
					//continue to next k-mer
					continue;
				}		
		}	
		
		reads[r].correct=false;
		break;
			
	}//if (bacs[0]!=bacSoFar) 
	//for (int p=0;p<POOLS;p++)
		//kmers[numKmers*POOLS+p]=e.getPoolFreq(p);
	//numKmers++;
}//end for k
}//end for r
//}//end while

for (int r=0;r<numReads;r++) 
	if (reads[r].correct)
		out<<reads[r].header<<endl<<reads[r].read<<endl;
	else 	
		out1<<reads[r].header<<endl<<reads[r].read<<endl;

out.close();
out1.close();
};

}

int main(int argc, char* argv[]) {

if (argc!=7)
	{
	std::cerr<<"\n\nUSE: hashfilter <path/ReadFile> <path/MapPoolsToBacsFile> <path/OutputFile> <k-mer_size> <input_file_number> <thread>\n\n";
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
for (unsigned int i=0;i<files;i++) {
	
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
