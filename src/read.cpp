#ifndef __EIGEN_H__
	#define __EIGEN_H__
	#include "constEigen.h"
#endif

#define PIX_SORT(a,b) { if ((a)>(b)) PIX_SWAP((a),(b)); }
#define PIX_SWAP(a,b) {float temp=(a);(a)=(b);(b)=temp; }

int readFloatFromFileEigen(FILE* fp, SpMatCol& phi, int m, int n) {
	int i,j;
	char* pch;
	char line[10000];
	
	std::vector<TripletF> tripletList;	
	
	for (i=0;i<m;i++) {
		if (fgets(line, 10000, fp) != NULL) { 
			pch = strtok (line," ");
			tripletList.push_back(TripletF(i,0,atof(pch)));
			for (j=1;j<n;j++) {
				pch = strtok (NULL, " \n");
				tripletList.push_back(TripletF(i,j,atof(pch)));
			}
		}	
	}
	phi.setFromTriplets(tripletList.begin(), tripletList.end());
	return 1;
};


void readBacMappings(unsigned short* bacPools) {
	int i, j, num;
	char* pch; 
	char line[1000];
	FILE* fp =fopen("/home/dduma/L1InfProjection/rice_mappings_sorted.txt", "r");
	if (fp == NULL) {
		fprintf(stderr, "Can't open mappings file! \n");
		exit(3);
	}
	for (i=0;i<BACS;i++) {
		//printf("%d\n", i);
		
		if (fgets(line, 1000, fp) != NULL) { 
			num=strlen(line);
			if (line[num-1]=='\n') {
				line[num-1]='\0';
				num--;
			}
			//printf("%s\n", line);
		}	
		
		//read the pool number
		pch=strtok(line, "\t");
		pch=strtok(NULL, "\t");
		pch=strtok(NULL, "\t");
		char* sig = pch;
		
		pch=strtok(sig, ",");
		bacPools[i*LAYERS + 0]=atoi(pch);
		//read in all 7 pools of this  bac
		for (j=1; j<LAYERS; j++) {
			pch=strtok(NULL,",");
			bacPools[i*LAYERS + j]=atoi(pch);
		}

	}//end for
	fclose(fp);

	//DEBUG
	/*for (i=0;i<BACS;i++) {
		printf("%d:\t", (i+1));
		for (j=0;j<LAYERS;j++)
			printf("%d ", bacPools[i*LAYERS + j]);
	printf("\n");
	}*/		
};


void readPoolMappings(int* poolBacs) {
	int j, num, pool;
	char *pch, *last; 
	
	char line[1000];
	FILE* fp =fopen("/home/dduma/L1InfProjection/pools2bacs.txt", "r");
	if (fp == NULL) {
		fprintf(stderr, "Can't open input file! \n");
		exit(3);
	}
	
	while (!feof(fp)) {
		if (fgets(line, 1000, fp)==NULL) 
			break;
		
		num=strlen(line);

		if (line[num-1]=='\n') {
			line[num-1]='\0';
			num--;
		}
		
		//read the pool number
		pch=strtok_r(line, ": ",&last);
		pool=atoi(pch);
		//if (pool==1)
		//	printf("%d: ", pool);

		//read in all the 169 bacs of pool
		for (j=0; j<BACS_IN_POOL; j++) {
			pch=strtok_r(NULL," ",&last);
			poolBacs[(pool-1)*BACS_IN_POOL + j]=atoi(pch);
			//if (pool==1)
			//	printf("%s,%d ", pch, poolBacs[j*m + (pool-1)]);

		}

	}//end while
	fclose(fp);

	//DEBUG
	/*for (i=0;i<m;i++) {
		printf("%d: ", (i+1));
		for (j=0;j<BACS_IN_POOL;j++)
			printf("%d ", poolBacs[j*m + i]);
	printf("\n");
	}*/		
};

int compareF(const void *_a, const void *_b) {  
  float *a, *b;      
  a = (float *) _a;
  b = (float *) _b;              
  //if(*a==*b)
  //  return 0;
  if(*a < *b)
    return -1;
  if(*a > *b)
    return 1;
  return 0;  
};

int compareSFRev(const void *_a, const void *_b) {
  struct residualF *a, *b;      
  a = (struct residualF *) _a;
  b = (struct residualF *) _b;
  
  //if((*a).norm==(*b).norm)
  // return 0;
  if((*a).norm > (*b).norm)
    return -1;
  if((*a).norm < (*b).norm)
    return 1;
  return 0;  
};

float findMedian(float a0, float a1, float a2, float a3, float a4, float a5, float a6) {

    PIX_SORT(a0, a5) ; PIX_SORT(a0, a3) ; PIX_SORT(a1, a6) ;
    PIX_SORT(a2, a4) ; PIX_SORT(a0, a1) ; PIX_SORT(a3, a5) ;
    PIX_SORT(a2, a6) ; PIX_SORT(a2, a3) ; PIX_SORT(a3, a6) ;
    PIX_SORT(a4, a5) ; PIX_SORT(a1, a4) ; PIX_SORT(a1, a3) ;
    PIX_SORT(a3, a4) ; return (a3) ;
};
