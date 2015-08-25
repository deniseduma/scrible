#include "classIO.h"

using namespace std;
using namespace general;
using namespace IOFn;

unsigned int hashSizes[SIZES]={5000549, 10000643, 20000003, 30000001, 40000003, 50000017, 60000011,70000027, 80000023, 90000049, 100000567, 150000001, 200000477, 250000013, 300000641, 350000041, 400000651, 450000007, 500000461, 550000001, 600000607, 650000011, 700000747, 750000007, 800000533, 850000013, 900000637, 950000017, 1000000531, 1250000027, 1500000713, 1750000027, 2000000579, 2147483647}; 

int main(int argc, char* argv[]) {

if (argc!=7)
	{
	std::cerr<<"\n\nUSE: tp <path/ReadFile> <path/MapPoolsToBacsFile> <path/OutputFile> <k-mer_size> <input_file_number> <thread>\n\n";
	exit(EXIT_FAILURE);
	}

unsigned short kmerSize=atoi(argv[4]);
unsigned short files=atoi(argv[5]); 
IOF<mybitsetx>* ioh=new IOF<mybitsetx>(argv[1], argv[2], argv[3], kmerSize, files);
ioh->ReadB();

ioh->computeProb();

return 0;
}
