#include "correction.h"

using namespace correction;

int main(int argc, char* argv[]) {

if (argc!=7)
	{
	std::cerr<<"\n\nUSE: corr <path/ReadFile> <path/MapPoolsToBacsFile> <path/OutputFile> <k-mer_size> <input_file_number> <thread>\n\n";
	exit(EXIT_FAILURE);
	}

string readFname(argv[1]);
unsigned int kmerSize=atoi(argv[4]);
unsigned int cores=atoi(argv[6]); 
unsigned int files=atoi(argv[5]); 

Timer* pTimer;

IOF ioh(argv[1], argv[2], argv[3], kmerSize, files);
ioh.ReadB();
//ioh.CheckCorrect();

//load BAC sigs into trie
unsigned short* bacPools=(unsigned short*)malloc(sizeof(unsigned short)*BACS*LAYERS);
readBacMappings(bacPools);

Trie trie;
for (unsigned short i=0;i<BACS;i++)
	trie.insert(&bacPools[i*LAYERS], 7, i+1);

HashTable<mysig> pools2bacs;
pools2bacs.allocHash();

for (unsigned int i=0;i<91;i++) {
	
	struct readFasta* reads=(struct readFasta*)malloc(MAX_POOL_SIZE*sizeof(struct readFasta));
	if (!reads) {
		perror("Error allocating the reads!");
		exit(1);
	}
	
	unsigned int numReads = ioh.readReads(i, reads);
	printf("Number of reads read is %d\n", numReads);  
	
	//correct reads in this pool
	MyCorrection* corr = new MyCorrection(&ioh, kmerSize, &trie, &pools2bacs, readFname, i, reads, numReads,cores);
	pTimer=new Timer("correct reads in pool");
	corr->correctReads();
	delete pTimer;
	delete corr;
	
	for (unsigned int r=0; r<numReads; r++) {
		free(reads[r].header);
		free(reads[r].read);
	}
	free(reads);

}//end for i

return 0;

}
