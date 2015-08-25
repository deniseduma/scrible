#include <omp.h>
#include <random>
#include <functional>
#include "classIO.h"

unsigned char c[8]={ 0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80};
unsigned char r[8]={ 0xFE, 0xFD, 0xFB, 0xF7, 0xEF, 0xDF, 0xBF, 0x7F};

unsigned int hashSizes[SIZES]={1000393, 2500009, 5000549, 10000643, 20000003, 30000001, 40000003, 50000017, 60000011,70000027, 80000023, 90000049, 100000567, 150000001, 200000477, 250000013, 300000641, 350000041, 400000651, 450000007, 500000461, 550000001, 600000607, 650000011, 700000747, 750000007, 800000533, 850000013, 900000637, 950000017, 1000000531, 1250000027, 1500000713, 1750000027, 2000000579, 2147483647}; 

using namespace std;
using namespace IOFn;

class IOF* p_ref;

//namespace IOFn {

//using namespace MYBIT;
//using namespace Eigen;

void showErrors(char* header, int* kmers, int numKmers, int pos, int p1, int p2) {

int nonZeros;
int numPools=0;

//count number of non-zero pools
for (int p=0;p<POOLS;p++) {
	nonZeros=0;
	for (int k=p1;k<p2;k++) 
		if (kmers[k*POOLS+p]>0) 
			nonZeros++;
	if (nonZeros>5)
		numPools++;
}

if (numPools>21)
	return; 

//print kmer number
printf("%s kmers %d errors pos=%d,p1=%d,p2=%d\n", header, numKmers, pos, p1, p2);
printf("   ");
for (int k=p1;k<p2;k++) 
	printf("%d ", k);
printf("\n");

for (int p=0;p<POOLS;p++) {
	//count number of non-zeros
	nonZeros=0;
	for (int k=p1;k<p2;k++) 
		if (kmers[k*POOLS+p]>0) 
			nonZeros++;
	
	if (nonZeros>5) {
		//print pool number
		if (p<10)
			printf(" %d: ", p);
		else
			printf("%d: ", p);
		
		//print kmers
		for (int k=p1;k<p2;k++) 
			printf("%d  ", kmers[k*POOLS+p]);
		printf("\n");
	}	
}		
printf("\n");	

};
//}

namespace IOFn {

void IOF::FastaToKmers(int readNo, char* header, char* read, int pos, int p1, int p2)
{

int num = strlen(read);
printf("%s\n", read);
int kmers[POOLS*MAX_KMERS];

int numKmers=0;
class mybitsetx b,d,e;
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
		
	if (!findN) {
		if (d<b) 
			e=d;
		else
			e=b;
		if (sh.searchCopy(e)) {
			for (int pool=0;pool<POOLS;pool++)
				kmers[numKmers*POOLS+pool]=e.getPoolFreq(pool);
			numKmers++;
		}
	}	
}//end for k

//deconv the overlap
showErrors(header, kmers, numKmers, pos, p1, p2);
};
}

void genErrors(int readNo, struct readFasta* read, 
	function <int ()> unisampler1, function <int ()> unisampler2) {
	int p1, p2;
	char r[100+1];
	strncpy(r, read->read, strlen(read->read));
	
	int pos = unisampler1();
	if (pos<26)
		p1=0;
	else if (pos>=26)
		p1=pos-25;
	p2=pos+1;
	
	p_ref->FastaToKmers(readNo, read->header, r, pos, p1, p2);
	char base=r[pos];
	int baseNo = (base=='A') ? 0 : (base=='C' ? 1 : (base=='G' ? 2 : 3));
	int newNo = unisampler2();
	while (newNo==baseNo)
		newNo=unisampler2();
	int newBase = (newNo==0) ? 'A' : (newNo==1 ? 'C' : (newNo=='2' ? 'G' : 'T'));
	r[pos]=newBase;
	p_ref->FastaToKmers(readNo, read->header, r, pos, p1, p2);
}


using namespace IOFn;

int main(int argc, char* argv[]) {

if (argc!=7)
	{
	std::cerr<<"\n\nUSE: hashfilter <path/ReadFile> <path/MapPoolsToBacsFile> <path/OutputFile> <k-mer_size> <input_file_number> <thread>\n\n";
	exit(EXIT_FAILURE);
	}

unsigned int qgram=atoi(argv[4]);
unsigned int cores=atoi(argv[6]); 
unsigned int files=atoi(argv[5]); 

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
	
	//generate random positions within the read
	uniform_int_distribution<int> uniform074(0,74);
	uniform_int_distribution<int> uniform03(0,3);
	mt19937 mtengine;
	//bind the RNG with the distribution
	auto unisampler1 = bind(uniform074,mtengine);
	auto unisampler2 = bind(uniform03,mtengine);

	int p1, p2, pos;
	#pragma omp parallel private(id,r) shared(reads,numReads)
	{
	id = omp_get_thread_num();
	#pragma omp for schedule(dynamic)
	for (r=0;r<10;r++) { 
		char read[100+1];
		strncpy(read, reads[r].read, strlen(reads[r].read));
		int pos = unisampler1();
		if (pos<26)
			p1=0;
		else if (pos>=26)
			p1=pos-25;
		p2=pos+1;
		
		p_ref->FastaToKmers(r, reads[r].header, read, pos, p1, p2);
		char base=read[pos];
		int baseNo = (base=='A') ? 0 : (base=='C' ? 1 : (base=='G' ? 2 : 3));
		int newNo = unisampler2();
		while (newNo==baseNo)
			newNo=unisampler2();
		char newBase = (newNo==0) ? 'A':(newNo==1 ? 'C':(newNo=='2' ? 'G':'T'));
		read[pos]=newBase;
		printf("old base %c, new base %c\n", base, newBase);
		p_ref->FastaToKmers(r, reads[r].header, read, pos, p1, p2);
		//genErrors(r, &reads[r], unisampler1, unisampler2);
	}	
	}
	
	for (r=0; r<numReads; r++) {
		free(reads[r].header);
		free(reads[r].read);
	}
	free(reads);

}//end for i

return 0;

}
