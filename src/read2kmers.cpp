#include<unistd.h> 
#include <stdio.h>

#ifndef _CIO_H_
	#define _CIO_H_
	#include "classIO.h"
#endif

namespace IOFn {

using namespace general;
using namespace std;
using namespace MYBIT;

void readData(vector<char*>& reads, FILE* in); 

template <class A_Type>
void IOFn::IOF<A_Type>::ReadToKmers(int id, int l,int u) {

//using namespace general;
//using namespace std;
//using namespace MYBIT;

FILE *in, *out;

clock_t startGlobal,endGlobal;

char delimC[] = ".";
Parser parser;
parser.update(delimC,ReadFname);
/*string fileExt="";
if (parser.size()>1)
	{
	fileExt=parser.get(1);
	}*/
for (int i=l;i<=u;i++)
{

startGlobal=clock();

ostringstream of,of1;
of<<parser.get(0)<<i<<"."<<parser.get(1);//<<"."<<parser.get(2)<<"."<<parser.get(3);
in=fopen(of.str().c_str(), "r");
if(in==NULL) 
	{
	cerr << "\nError opening input file " << of.str() << endl;
	exit(EXIT_FAILURE);
	}

vector<char*> reads;
readData(reads, in);
fclose(in);

of1<<OFname<<i<<".kmers" ;
out=fopen(of1.str().c_str(),"w");
if (out==NULL)
	{
	cerr << "\nError opening ouput file " << of1.str() << endl;
	exit(EXIT_FAILURE);
	}
cout << "[ReadToKmers]Thread " << id << " pool " << i << endl;
cout << "[ReadToKmers]Input file is " << of.str() <<  "\nOutput file is " << of1.str() << endl << endl;


//omp_set_num_threads(1);
//#pragma omp parallel shared(reads, out)  
//{
//printf("thread %d numReads %lu\n", omp_get_thread_num(), reads.size());
//#pragma omp for schedule(dynamic)
for (unsigned int r=0; r<reads.size(); r+=2)  {
	
	char* name=reads[r]; 
	char* read=reads[r+1];
	unsigned int num=strlen(read);
	
	if ((num>=kmerSize+1)&&((read[0]=='A')||(read[0]=='G')||(read[0]=='C')||(read[0]=='T')||(read[0]=='N'))) {
	
	vector<vector<unsigned short> > kmers;
	vector<unsigned short> kmersInPool(POOLS, 0);

	class mybitsetx b,d,e;
	for(unsigned int k=0; k<=num-(kmerSize)+1; k++) {
		int z=0;
		int j=DIM*2-1;
		bool findN=false;
		findN=false;
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
		}		
		
		if (!findN) {
			if (d<b) 
				e=d;
			else
				e=b;
			
			if (sh.searchCopy(e)) {
				vector<unsigned short> kmer(POOLS);
				for (int pool = 0; pool<POOLS; pool++) {
					kmer[pool]=e.getPoolFreq(pool);
					//kmer[pool]=e.getPool(pool);
					kmersInPool[pool]+=e.getPool(pool);
				}
				kmers.push_back(kmer);
			}
		}
	}//end for k
	
	fprintf(out, "%s kmers: %lu\n", name, kmers.size());
	for (unsigned int pool=0; pool<POOLS; pool++) {
		if (kmersInPool[pool]>0) {
			fprintf(out, "%d ", pool);
			for (unsigned int k=0; k<kmers.size(); k++) {
				fprintf(out, "%d ", kmers[k][pool]);
			}	
			fprintf(out, "\n");
		}
	}	
	
	}//end if 

}//end for r
//}//end pragma omp

for (unsigned int j=0; j<reads.size(); j++) {
	delete [] reads.at(j);
}

fclose(out);
endGlobal=clock();

cout<<"\n\tTime to output the kmers in pool "<<i<<": "<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl<<endl;
		
}//end for(i=l;i<u;i++)

};

void readData(vector<char*>& reads, FILE* in) 
{
int num;
char* buffer;	
while (!feof(in)) {
	buffer=new char[MAXSIZE];	
	if (fgets(buffer, MAXSIZE, in)==NULL) {
		break;
	}	
	num=strlen(buffer);
	if (buffer[num-1]=='\n')
	{
		buffer[num-1]='\0';
		num--;
	}
	reads.push_back(buffer);
}
	//DEBUG
	printf("[readData]numReads read from file %lu\n", reads.size());
};


/*Builds hashtable on all reads (after correcting them using sga)
and resizes the hashtable if necessary. 
Avoids ReadLight() followed by Read() which is extremely time consuming.*/

/*template <class A_Type>
void IOFn::IOF<A_Type>::ReadBarley(){

//using namespace general;
//using namespace std;
//using namespace MYBIT;

clock_t startGlobal,endGlobal;

unsigned int num;
char buffer[MAXSIZE];

//ifstream in;
FILE* in;

char delimC[] = ".";
Parser parser;
parser.update(delimC,ReadFname);

for (unsigned int i=0;i<numFiles;i++)
{
//DEBUG
cout << endl << "***Pool " << i << "***" << endl;

ostringstream of;
of<<parser.get(0)<<i<<"."<<parser.get(1);//<<"."<<parser.get(2)<<"."<<parser.get(3);

//in.open(of.str().c_str(),ofstream::in);
in=fopen(of.str().c_str(), "r");
if(in==NULL) 
	{
	cerr << "\nError opening input file " << of.str() << endl;
	exit(EXIT_FAILURE);
	}

startGlobal=clock();
while (!feof(in))
	{
	buffer[0]='\0';
	if (fgets(buffer, MAXSIZE, in)==NULL)
		break;
	//if (ferror(in))
	//	perror("the bug is");
	//num=in.gcount();
	num=strlen(buffer);
	*if (buffer[num-1]!='\0')
		{
		buffer[num]='\0';
		num++;
		}*
	if (buffer[num-1]=='\n') {
		buffer[num-1]='\0';
		num--;
	}
	if ((num>=kmerSize+1)&&((buffer[0]=='A')||(buffer[0]=='G')||(buffer[0]=='C')||(buffer[0]=='T')||(buffer[0]=='N')))
		{
			readCounts[i]++;	
			insert(buffer,num,i);
			*if(numReads%20000000==0) {
				cout<<"[Read]reads so far: "<<numReads<<" : "<<" stored kmers: "<<sh.size()<<" rejected  kmers: "<<numDiscard<<endl;
				//sh.print();
				sh.info();
			}	*
		}
	}//end while (!in.eof())

numReads+=readCounts[i];

fclose(in);

endGlobal=clock();
cout << "\n[classIO:Read]Time to read pool "<<i<<": "<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl;
cout<<"[classIO:Read]Reads in this pool: "<<readCounts[i]<<" kmers: "<<kmerCounts[i]<<" kmers stored so far: "<<sh.size() << " k-mers discarded " << numDiscard <<endl;

};

//DEBUG
cout <<"\n[classIO:Read]Done reading " << numFiles << " pools for a total of " << numReads << " reads." << endl;
cout << "\n[classIO:Read]After reading all pools sh has size " << sh.size() << " and capacity " << sh.capacity() << endl << endl; 
sh.checkCorrect();

};*/

/*template <class A_Type>
int IOFn::IOF<A_Type>::readReads(unsigned int i, struct readFasta* reads) {

//using namespace general;
//using namespace std;
//using namespace MYBIT;

ostringstream of;

char delimC[] = ".";
Parser parser;
parser.update(delimC, ReadFname);

ifstream in;
of<<parser.get(0)<<i<<"."<<parser.get(1);
in.open(of.str().c_str(), ios_base::in);
if(!in) {
	cerr<<"Error opening input file "<<of.str()<<endl;
	exit(EXIT_FAILURE);
}

//DEBUG
cout << "Reading file " << of.str() << endl;

unsigned int num;
char head[MAXSIZE];
char seq[MAXSIZE];

int index=-1;
int numReads=0;
struct readFasta* crt;

clock_t startGlobal, endGlobal;

startGlobal=clock();
while (in.good()) {
	
	head[0]='\0';
	in.getline(head, MAXSIZE);
	num=strlen(head);
	
	if ((num>0)&&((head[0]=='@')||(head[0]=='>'))) {
		
		index++;
		crt = reads + index;
		
		numReads++;
	
		crt->correct=true;
		crt->header=(char*)malloc((num+1)*sizeof(char));
		strncpy(crt->header, head, num);
		crt->header[num]='\0';
	
		//DEBUG
		//printf("%s\n", crt->header);
		//printf("numOverlaps so far %d\n", numOverlaps);
	
		seq[0]='\0';
		in.getline(seq, MAXSIZE);
		num=strlen(seq);
					
		assert((num>=75)&&((seq[0]=='A')||(seq[0]=='G')||(seq[0]=='C')||(seq[0]=='T')||(seq[0]=='N'))); 
		
		crt->read=(char*)malloc((num+1)*sizeof(char));
		strncpy(crt->read, seq, num);
		crt->read[num]='\0';
			
	}//end if
}//end while

in.close();

endGlobal=clock();
cout << "\n[read2kmers:readReads]Time to read pool "<<i<<": "<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl;
	
return numReads;

};*/


/*void IOF::FastaToKmers(int readNo, char* header, char* read, int pos, int p1, int p2)
{

using namespace Eigen;

//char* header=r->header;
//char* read=r->read;
int num = strlen(read);

//SpMatCol kmers(POOLS, MAX_KMERS);
//vector<TripletF> tripletList(POOLS*MAX_KMERS);
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
				//tripletList.push_back(TripletF(pool,numKmers,e.getPoolFreq(pool)));
				kmers[numKmers*POOLS+pool]=e.getPoolFreq(pool);
			numKmers++;
		}
	}	
}//end for k
	
//kmers.setFromTriplets(tripletList.begin(), tripletList.end());
printf("numKmers=%d\n", numKmers);
//deconv the overlap
//showErrors(header, kmers, numKmers, p1, p2);
};*/

}
