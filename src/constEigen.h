#include <Sparse>
#include <Dense>
#include <set>

#ifndef __CON_H__
	#define __CON_H__
	#include "const.h"
#endif

//#define POOLS2C 12
//#define POOLS 91
//#define BACS 2197
//#define MAX_BACS 30
//#define LAYERS 7
//#define BACS_IN_POOL 169
//#define POOLS_PER_LAYER 13

//#define MAX_POOL_SIZE 18000000 //for barley

#define CORES 10

#define MAX_DIFF_X 1e-6
#define MAX_DIFF_OBJ 1e-2
#define MIN_ITER 10

//0 direct, 1 iterative
#define METHOD 3
#define MAX_TOL 1e-6
#define MAX_ITER 91

#define T 1
//#define TOP 3
#define SPARSITY 3
#define KMER_TRESHOLD 2
#define BAC_SCORE_TRESHOLD 0.5
#define DIST_FROM_MAX 0.5
#define KEEP_IN_COL 21
#define KEEP_IN_LAYER 3 
#define LAYERS_2_MATCH 7 
#define RESIDUAL_FROBENIUS_NORM 1e-6
#define SCORE_FUNC 0 

typedef Eigen::SparseMatrix<float, Eigen::ColMajor> SpMatCol; 
typedef Eigen::SparseMatrix<float, Eigen::RowMajor> SpMatRow; 
typedef Eigen::Triplet<float> TripletF;

struct residualF{
    int rowInd;
    float norm0;
    float norm;
};

struct residualD{
    int rowInd;
    double norm0;
    double norm;
};

struct layer_set{
	int count;
	//int pools[POOLS_PER_LAYER];
	//float vals[POOLS_PER_LAYER];
	int pools[KEEP_IN_COL];
	float vals[KEEP_IN_COL];
};

struct readKmers {
	//public: 
	char* header;
	int numKmers;
	int numBacs;
	int* bacs; 
	SpMatCol* pools;
};

struct readFasta {
	char* header;
	char* read;
	bool correct;

	readFasta(): header(new char[MAXSIZE]),read(new char[MAXSIZE]), correct(false){}

	~readFasta() {delete [] header; delete [] read;}
};

struct overlap {
	char* header;
	int numReads;
	char** reads;
};

struct kmerInfo {
	unsigned int pos;
	unsigned int numPools;
	unsigned short pools[POOLS];
	unsigned int numBacs;
	unsigned short bacs[MAX_BACS];
};

/*struct thread_data{
	int l;
	int u;
};*/

struct cset {
	//bool disagree;
	short begin, end;
	//short initBegin, initEnd;
	std::set<unsigned short> bacs;
	double centroid[POOLS];
	unsigned short numPools;
		
	cset(): begin(-2), end(-2), numPools(0) {
		//for (unsigned short i=0;i<POOLS;i++)
		//	centroid[i]=0;
	};
	
	cset(const struct cset& b): begin(b.begin), end(b.end), bacs(b.bacs), numPools(b.numPools) {
		for (unsigned short i=0;i<POOLS;i++)
			centroid[i]=b.centroid[i];
	};

	struct cset& operator=(const struct cset& b) {
		if (this!=&b) {
			begin=b.begin;
			end=b.end;
			bacs=b.bacs;
			for (unsigned short i=0;i<POOLS;i++)
				centroid[i]=b.centroid[i];
			numPools=b.numPools;	
		}	
		return *this;
	};
	
	bool operator==(const struct cset& b) {
		if ((begin==b.begin)&&(end==b.end))
			return true;
		return false;	
	};

	/*void initEnds() {
		initBegin=begin;
		initEnd=end;
	};*/
};


double normL1Inf(double *m, int nRows, int nCols);
float normFrobenius(float *m, int nRows, int nCols);
int compareF(const void *_a, const void *_b); 
int compareD(const void *_a, const void *_b); 
int compareFRev(const void *_a, const void *_b);
int compareDRev(const void *_a, const void *_b);
int compareSF(const void *_a, const void *_b);
int compareSD(const void *_a, const void *_b);
int compareSFRev(const void *_a, const void *_b);
int compareSDRev(const void *_a, const void *_b);

/*void transpose(int* A, int* ATrans, int m, int n);
void inverse(int* A, int* AInv, int n);
void matrixMultiply(int* A, int* B, int* C, int m, int n, int p);

void projL1Inf(double B[], double C, double A[], double w[], int nRows, int nCols);
double computeError(int* A, int* M, int* X, int* G,  int m, int n, int q, struct residualD* maxes, int l);

void buildMatrixBis(struct read* r, double** pM, int m, int q); 
int readIntFromFile(FILE* fp, int**  A, int m, int n);
int readIntFromFileRowMajor(FILE* fp, int**  A, int m, int n);
*/
//int readFloatFromFile(FILE* fp, float**  A, int m, int n);
int readFloatFromFileEigen(FILE* fp, SpMatCol&  phi, int m, int n);
//int readData(struct read* reads, FILE* fp);
//int readDataWithBacs(struct read* reads, FILE* fp);

void readBacMappings(unsigned short* bacPools);
void readPoolMappings(int* poolBacs);

/*void runOpt(int* A, int*ATrans, int* M, double *X, double* G, double* R,
	double* w, int m, int n, int q, double eta0, double C); 

double computeError(int* A, int* M, int* X, int* G,  int m, int n, int q, struct residualD* maxes, int l);

void printMatrix(double *m, int nRows,int nCols);
void printMatrixInt(int* A, int nRows,int nCols);
void printMatrixDouble(double* A, int nRows,int nCols);
void printMatrixFloat(float* A, int nRows,int nCols);
void printResult(int* A, int* M, double* X, int m, int n, int q, char* read, int* bacs, FILE* out);	

float quick_select(float arr[], int n); 
float findMedian(float a0, float a1, float a2, float a3, float a4, float a5, float a6);
*/

/*class myhash
{
public:
	size_t operator()(const mybitsetx &b) const {
		return b.trans2();
	}
};


class myequal_to
{
public:
	bool operator()(const mybitsetx& b1, const mybitsetx& b2) const {
		return (b1==b2);
	}
};

class mycompareRev
{
public:
	bool operator()(unsigned int c1, unsigned int c2) const {
		if (c1>c2) 
			return true;
		return false;	
	}
};

typedef multimap<unsigned int, mybitsetx, mycompareRev> my_m_map;
typedef unordered_map<mybitsetx, unsigned int, myhash, myequal_to> my_u_map;
}*/
