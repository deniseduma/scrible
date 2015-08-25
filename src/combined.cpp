#include <stdio.h>
#include <unistd.h> 
#include <math.h>
#include <vector>

#include <assert.h>

//#include "constEigen.h"
//#include <LU>
//#include <Cholesky>
//#include <Eigenvalues>
//#include <SparseCholesky>
//#include <IterativeLinearSolvers>
#include "classIO.h"

using namespace std;
using namespace MYBIT;
using namespace general;
using namespace IOFn;
using namespace Eigen;


unsigned char c[8]={ 0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80};
unsigned char r[8]={ 0xFE, 0xFD, 0xFB, 0xF7, 0xEF, 0xDF, 0xBF, 0x7F};

unsigned int hashSizes[SIZES]={1000393, 2500009, 5000549, 10000643, 20000003, 30000001, 40000003, 50000017, 60000011,70000027, 80000023, 90000049, 100000567, 150000001, 200000477, 250000013, 300000641, 350000041, 400000651, 450000007, 500000461, 550000001, 600000607, 650000011, 700000747, 750000007, 800000533, 850000013, 900000637, 950000017, 1000000531, 1250000027, 1500000713, 1750000027, 2000000579, 2147483647}; 

/*namespace IOFn {

void IOF::deconvOverlap(const SpMatCol& phi,int* bacPools,char* header,int numKmers,const SpMatCol& Y, ofstream& out) {

int t, i, j, s;

//int numBacs=read->numBacs;
//int* bacs=read->bacs;
int numBacs=BACS;
int bacs[numBacs];
for (i=0;i<numBacs;i++)
	bacs[i]=i+1;

//DEBUG
//printf("header %s\n", header);
//printf("numKmers %d\n", numKmers);
*printf("matrix Y\n");
for (i=0;i<POOLS;i++) {
	printf("%d: ", i); 
	for (j=0;j<numKmers;j++)
		printf("%.0f ", Y.coeff(i,j));
	printf("\n");
}*

SpMatCol R(Y);

float median;
//float matched[LAYERS];

int supportSize=0;
int support[numBacs];
struct residualF scores[numBacs];

struct layer_set top_pools[LAYERS];
struct residualF maxes[POOLS];
	
std::vector<TripletF> tripletList;
//tripletList.reserve(POOLS*T*SPARSITY);
//SpMatCol ARed;
//ARed.resize(POOLS, T*SPARSITY);
//SpMatCol Z(BACS, numKmers);

float tol, resid;
float normY = Y.norm();
if (normY == 0)
	normY = 1;

//DEBUG
//printf("t=%d, normR=%f, normW=%f\n", 0, R.norm(), normW);

//while (R.norm()>0.25*normW) {
for (t=1;t<=T;t++) {
	
	if (numKmers<3)
		break;
		
	int iter = MAX_ITER;
	
	//initialize bac scores
	for (int b=0;b<numBacs;b++) {
		scores[b].rowInd=b;
		scores[b].norm0=0;
		scores[b].norm=0;
	}
	
	for (int k=0; k<R.outerSize(); k++) {//for each k-mer
	
		int numPools=0;
		//keep top x pool measurements in each layer
		for (SpMatCol::InnerIterator it(R, k); it; ++it) {
			maxes[numPools].rowInd=it.row();
			maxes[numPools++].norm=it.value();
		}
		qsort(maxes, numPools, sizeof(struct residualF), &compareSFRev);
	
		for (j=0;j<LAYERS;j++) {
			top_pools[j].count=0;
		}
		
		for (j=0;j<numPools;j++) {
			int pool=maxes[j].rowInd;
			float val = maxes[j].norm;
			int layer=pool/POOLS_PER_LAYER;
			if (val>0 && top_pools[layer].count<KEEP_IN_LAYER) {
				top_pools[layer].pools[top_pools[layer].count]=pool;
				top_pools[layer].vals[top_pools[layer].count]=val;
				top_pools[layer].count++;
			}
		}//end for j
	
		//assert(numBacs>3);

		for (int b=0;b<numBacs;b++) {//for each bac 
			if (bacs[b]==0)
				continue;
			
			int bac=bacs[b] - 1;
			int* pools=&bacPools[bac*LAYERS];
			
			int layers_matched=0;
			for (j=0;j<LAYERS;j++) {//for each layer
				int pool=pools[j]-1;
				for (i=0;i<top_pools[j].count;i++) {//for each retained pool in this layer
					if (pool==top_pools[j].pools[i]) {
						//matched[layers_matched]=top_pools[j].vals[i];
						layers_matched++;
						break;
					}
				}//end for i	
			}//end for j
		
			if (layers_matched>=LAYERS_2_MATCH) {
				//using median instead
				//the quickselect method
				//median=quick_select(matched, layers_matched);
				
				//update the bac score
				if (SCORE_FUNC==0) {
					//if (fabs(min)>=KMER_TRESHOLD) {
					scores[b].norm0++;
					scores[b].norm++;
					//}
				} else if (SCORE_FUNC==1) {	
					//if (fabs(median)>=KMER_TRESHOLD) {
					scores[b].norm0++;
					scores[b].norm += fabs(median);
					//}
				} else if (SCORE_FUNC==2) {
					//if (fabs(median)>=KMER_TRESHOLD) {
					scores[b].norm0++;
					scores[b].norm += median*median;
					//}
				}
			}//end if layers matched
		}//end for b
	}//end for k
	
	//sort bacs in reverse order of their scores
	qsort(scores, numBacs, sizeof(struct residualF), &compareSFRev);
	
	//DEBUG
	//for (s=0;s<5;s++)
	//	printf("%d %.0f\n", bacs[scores[s].rowInd] - 1, scores[s].norm);
	
	//return the top l bacs according to the score function
	for (s=0;s<SPARSITY;s++) {
		if (scores[s].norm0<BAC_SCORE_TRESHOLD*numKmers) 
			break;
		if (scores[s].norm0<DIST_FROM_MAX*scores[0].norm0)  
			break;
		support[supportSize++]=bacs[scores[s].rowInd] - 1;
		bacs[scores[s].rowInd]=0;
	}	
	
	if (s==0) 
		break;
	
	//project onto the span of atoms chosen so far
	
	SpMatCol ARed;
	if (METHOD==0 || METHOD==1 || METHOD==2) {
		ARed.resize(POOLS, supportSize);
	
		for (i=0;i<s;i++) {
			int col=support[supportSize-s+i];
			for (SpMatCol::InnerIterator it(phi,col); it; ++it) 
				tripletList.push_back(TripletF(it.row(),supportSize-s+i,it.value()));
		}
	
		ARed.setFromTriplets(tripletList.begin(), tripletList.end());
	
	}	
	
	if (METHOD==0) {//DIRECT METHOD
		SpMatCol AT=ARed.transpose();
		SpMatCol ATA=AT*ARed;
	
		SimplicialLDLT<SpMatCol> solver;
		//solver.compute(ARed);
		solver.compute(ATA);
		//SpMatCol Z(BACS, numKmers); 
		//Z = solver.solve(M);
		SpMatCol S(BACS, POOLS); 
		S = solver.solve(AT);
		SpMatCol P=ARed*S;
		//SpMatCol PW = P*W;
		R=Y-P*Y;
	} else if (METHOD==1) {//RICHARDSON ITERATION
	 
	 	//estimate of largest singular value  of ARed
	 	float coherenceOfA=2.0/7.0;
		float sigma = coherenceOfA * supportSize;
		
		SpMatCol Z;
		Z.resize(supportSize, numKmers);
		//SpMatCol ARZ = ARed*Z;
		SpMatCol Res = Y - ARed * Z;
								
		if ((resid = Res.norm() / normY) <= MAX_TOL) {
	        	tol = resid;
	        	iter = 0;
		 } else {
			//using the transpose of ARed
			SpMatCol AT = ARed.transpose();
			AT /= sigma*sigma; 

			for (int it = 0; it < MAX_ITER; it++) {
				//SpMatCol ATR = AT*Res;
				Z = Z + AT * Res;
								
				//iterate
				//SpMatCol ARZ = ARed*Z;
				Res = Y - ARed * Z;
				//DEBUG
				printf("\tresid norm %f, Y norm %f, resid %f, MAX_TOL %f\n\n", Res.norm(), normY, Res.norm()/normY, MAX_TOL);  
		        
				if ((resid = Res.norm() / normY) <= MAX_TOL) {
					tol = resid;
					iter = it;
					break;
			    	}
			}//end for it
		}//end else	
		//update residual
		//SpMatCol AZ=A*Z;
		R=Y - phi * Z;

		//DEBUG
		printf("END %s\n", header);
		printf("iter %d, tol %f\n\n", iter, tol);
	
	} else if (METHOD==2) {//the method of normal equations
		SpMatCol AT=ARed.transpose();
		SpMatCol ATA=AT*ARed;
		SpMatCol ATW=AT*Y;
		
		SpMatCol Z;
		Z.resize(supportSize, numKmers);
		
		//solve least square using conjugate gradient
		//SimplicialLDLT<SpMatCol> solver;
		ConjugateGradient<SpMatCol> solver;
		//BiCGSTAB<SpMatCol> solver;
		solver.compute(ATA);
		Z = solver.solve(ATW);
				
		SpMatRow Zf;
		Zf.resize(BACS, numKmers);
		SpMatRow Zr(Z);
		
		std::vector<TripletF> tList;
		tList.reserve(supportSize*POOLS);
		
		for (i=0;i<Zr.outerSize();i++) {
			for (SpMatRow::InnerIterator it(Zr, i);it;++it ) {
				//Zf.coeffRef(support[i], col)=Zr.coeff(i,col);
				tList.push_back(TripletF(support[i], it.col(), it.value()));
			}
		}
		Zf.setFromTriplets(tList.begin(), tList.end());
		//SpMatCol Zfc(Zf);
		R = Y - phi * SpMatCol(Zf);
		
	}//end METHOD 2
	
}//end for t	
//} end while (R.norm > 0.25*normW) 

int num;
char str[100];
strcat(header, " bacs: ");
for(i=0;i<supportSize;i++) {
	str[0]='\0';
	num=sprintf(str, "%d ", support[i] + 1);
	strncat(header, str, num);
}	
out<<header<<endl;
};
}
*/

int main (int argc, char** argv) {

if (argc!=7) {
	cerr<<"Usage: combined <path/ReadFile> <path/MapPoolsToBacsFile> <path/OutputFile> <k-mer_size> <input_file_number> <thread number>\n";
	exit(1);
}

char* inFile=argv[1];
char* outFile=argv[3];
int kmerSize=atoi(argv[4]); 
int numFiles=atoi(argv[5]);
int cores=atoi(argv[6]);

class IOF ref(inFile, argv[2], outFile, kmerSize, numFiles);
ref.ReadB();

time_t t0, t1;
clock_t c0, c1;

FILE* fp = fopen("/home/dduma/L1InfProjection/pools_91_2197Sorted.int", "r");
if (fp == NULL) {
	fprintf(stderr, "Can't open design file!\n");
	exit(1);
};

SpMatCol phi(POOLS, BACS);
//A.reserve(VectorXi::Constant(BACS,LAYERS));
t0=time(NULL);
c0=clock();
readFloatFromFileEigen(fp, phi, POOLS, BACS);
t1=time(NULL);
c1=clock();
fclose(fp);
printf("time to read matrix phi is %.2f %.2f\n", (float)(c1 - c0)/CLOCKS_PER_SEC, (float)(t1-t0));
//normalize A

phi /= 7;

int* bacPools=(int*)malloc(sizeof(int)*BACS*LAYERS);
readBacMappings(bacPools);

int id, r;
struct overlap* overlaps;

for (int i=0;i<numFiles;i++) {
	
	overlaps=(struct overlap*)malloc(MAX_POOL_SIZE*sizeof(struct overlap));
	if (!overlaps) {
		perror("Error allocating overlaps: ");
		exit(1);
	}
	
	ofstream out;
	int numOverlaps = ref.readOverlaps(i, overlaps, out);
	printf("Number of overlaps read is %d\n", numOverlaps);  
	
	omp_set_num_threads(cores);
	printf("Number of cores for deconvolution is %d\n", cores);	
	
	t0=time(NULL);
	c0=clock();
	#pragma omp parallel private(id,r) shared(phi,overlaps,numOverlaps,bacPools,out)
	{
	id = omp_get_thread_num();
	#pragma omp for schedule(dynamic)
	for (r=0;r<numOverlaps;r++) {
		ref.OverlapToDeconv(id, r, phi, bacPools, &overlaps[r], out);
	}
	}
	t1=time(NULL);
	c1=clock();
	printf("Time to deconvolute reads is %.2f %.2f\n",(float)(c1 - c0)/CLOCKS_PER_SEC, (float)(t1-t0));
	
	out.close();
	
	for (r=0; r<numOverlaps; r++) {
		free(overlaps[r].header);
		for (int s=0;s<overlaps[r].numReads;s++)
			free(overlaps[r].reads[s]);
		free(overlaps[r].reads);
	}
	free(overlaps);

}//end for i

free(bacPools);

return 1;

}
