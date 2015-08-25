#ifndef __EIGEN_H_
	#define __EIGEN_H__
	#include "constEigen.h"
#endif

#ifndef __FSTREAM__
	#define __FSTREAM__
	#include <fstream>
#endif

#ifndef __SSTREAM__
	#define __SSTREAM__
	#include <sstream>
#endif 
#ifndef __IOS_H__
	#define __IOS_H__
	#include <iostream>
#endif 

#ifndef __STDL__
	#define __STDL__
	#include <stdlib.h>
#endif

#ifndef __STR_H__
	#define __STR_H__
	#include <string.h>

#endif
 
#ifndef __TIM_H__
	#define __TIM_H__
	#include <time.h>

#endif

#include <set>
#include <map>
//#include <multimap.h>
#include <unordered_map>

#include <omp.h>

#ifndef HASH_H_
	#define HASH_H_
	#include "hashChaining.h"
#endif

#ifndef __GNR_H__
	#define __GNR_H__
	#include "general.h"

#endif

//#include <set>
#include <map>
//#include <multimap.h>
#include <unordered_map>

#define DEBUGBAC 0

extern unsigned int hashSizes[SIZES];

namespace IOFn
{
using namespace std;
using namespace MYBIT;
using namespace Cl_HASH;

template <class A_Type> 
class IOF
{
//!Input file, this name is concatented with an integer, so we assume the input file name [A-Ba-b0-9]+[0-1]+
string ReadFname;
//!Input file for mapping pools into a BAC.
string BACFname;
//!Output file
string OFname;
//!kmerSize size
unsigned int kmerSize;
//!number of input files
unsigned int numFiles;
//!Read number
unsigned int numReads;
//!Max number of reads in any pool
unsigned long long maxReads;
//!Max number of kmers in any pool
unsigned long long maxKmers;
//!Reject kmerSizes number
unsigned int numDiscard;
//!Per pool read counts
unsigned long long readCounts[POOLS];
//!Per pool kmer counts
unsigned long long kmerCounts[POOLS];
//!kmer hash table 
HashTable<A_Type> sh;
//! BAC hash table 
HashTable<mybitsetBAC> shBAC;
//!Hash table size
//short crtHS;

public:

//!Construtor. It takes in input/output file names, kmerSize size and the number of input files
IOF(const char* RFname,const char* BFname,const char* OFname, unsigned int kmerSize, unsigned int numFiles):
ReadFname(RFname), BACFname(BFname), OFname(OFname), kmerSize(kmerSize), numFiles(numFiles), numReads(0), maxReads(0), maxKmers(0), numDiscard(0) {
	for (unsigned int i=0; i<POOLS;i++) {
		readCounts[i]=0;
		kmerCounts[i]=0;
	}	
	//memset(readCounts,0, POOLS*sizeof(unsigned long long) );
	//memset(kmerCounts,0, POOLS*sizeof(unsigned long long) );
	cout << "Constructor of IOF called" << endl;
	cout << "ReadFname is " << ReadFname << endl;
	cout << "Total size so far is " << sizeof(sh) + sizeof(shBAC) + 2*POOLS*sizeof(unsigned long long) << endl;
};


//!Construtor. It takes in input/output file names.
IOF(const char* RFname,const char* BFname,const char* OFname): 
ReadFname(RFname), BACFname(BFname), OFname(OFname), kmerSize(DIM), numFiles(POOLS), numReads(0), maxReads(0), maxKmers(0), numDiscard(0) {
	for (int i=0; i<POOLS;i++) {
		readCounts[i]=0;
		kmerCounts[i]=0;
	}	
	//memset(readCounts,0, POOLS*sizeof(long long) );
	//memset(kmerCounts,0, POOLS*sizeof(long long) );
	cout << "Constructor of IOF called" << endl;
	cout << "ReadFname is " << ReadFname << endl;
	cout << "Total size so far is " << sizeof(sh) + sizeof(shBAC) + sizeof(readCounts) + sizeof(kmerCounts) << endl;
};


//!It takes in input the hash table size and allocates it.
void AllocHash(){
	cout << "[AllocHash]Allocating hashtables." << endl;
	//sh.allocHash(5);
	sh.allocHash(3);
	shBAC.allocHash(0);	
};

//!Destructor
~IOF() {
	//delete sh;
	//delete shBAC;
	cout << "Destructor of IOF called" << endl;
}

//!It reads the input files containing reads and stores the their data into the Hash table 
void Read();
void RetrieveBacs(Trie* trie) {
	cout << "[RetrieveBacs]Retrieving BAC signatures." << endl;
	sh.RetrieveBacs(trie);
}
void ReadBarley();
//!It reads the input files containing reads and stores the their data into the Hash table 
void ReadLight();
void ReadIntermediate();
//!It writes the hash table data in the output file
void Write();
//!It reads the input files cantaining the BACs and stores the their data into the Hash table 
void ReadBac();
//!It returns the size of the k-mer hash table
unsigned int size() {return sh.size();}
//!Compare new hash table with the old one to make sure everything is correct
void compareTables();
/*{
	using namespace general;
	Cl_HASH::compareTables(shOld, sh);
}*/	

//!It reads the input files (in binary format) storing their data into the hash table 
inline bool ReadB(){
	using namespace general;
	FILE * pFile;
	char delimC[] = ".";
	Parser parser;
	parser.update(delimC,ReadFname);
	string tmp=parser.get(0)+".count";
  	pFile = fopen ( tmp.c_str(), "r" );
	if (pFile==NULL)
  		{
		cerr<<"\n" << tmp << " not found!" << endl;
		cerr<<"Hash Table binary encoding not found!"<<endl <<endl;
		return false;
		}
	int POOLStmp=0;
	if (fread (&POOLStmp, sizeof(unsigned int),1, pFile )!=	1)	
		{
		cerr<<"Error reading POOLStmp"<<endl; 
		exit (EXIT_FAILURE);
		}
	if (POOLStmp!=POOLS)
		{
		cerr<<"\n*****Error you read a different number of pools ("<<POOLStmp<<") *****" << endl;
		return false;
		}
	if (fread (readCounts, sizeof(unsigned long long),POOLS, pFile )!=POOLS)
		{
		cerr<<"Error reading readCounts"<<endl; 
		exit (EXIT_FAILURE);
		}
	if (fread (kmerCounts, sizeof(unsigned long long),POOLS, pFile )!=POOLS)
		{
		cerr<<"Error reading kmerCounts"<<endl; 
		//perror("Error reading kmerCounts");
		exit (EXIT_FAILURE);
		}
	
	for (int i=0;i<POOLS;i++) {
		if (maxReads<readCounts[i]) 
			maxReads=readCounts[i];
		if (maxKmers<kmerCounts[i])
			maxKmers=kmerCounts[i];
	}

	cout << "[ReadB]maxKmers is " << maxKmers << endl;
	
	bool out1=sh.Read(parser.get(0));
	return ((out1)) && (shBAC.Read(BACFname));

};

//!It writes the hash table data in the output file (in binary format).
inline void WriteB(){
	using namespace general;
	FILE * pFile;
	char delimC[] = ".";
	Parser parser;
	parser.update(delimC,ReadFname);
	string tmp=parser.get(0)+".count";
  	pFile = fopen ( tmp.c_str(), "w" );
	if (pFile==NULL)
  		{
		cerr<<"\nError opening output file " << tmp  << endl;
		exit(EXIT_FAILURE);
		}
	//write pool number
	int POOLStmp=POOLS;
	fwrite (&POOLStmp, sizeof(int), 1, pFile );
	//write number of reads per pools
	fwrite (readCounts, sizeof(unsigned long long), POOLS, pFile );
	//write number of kmers per pools
	fwrite (kmerCounts, sizeof(unsigned long long), POOLS, pFile );
	fclose (pFile);
	cout<<"[classIO:WriteB]Per pool read and kmer frequencies saved in: "<<tmp<<endl;
	sh.Write(parser.get(0));
	shBAC.Write(BACFname);

};

void  print(){
	sh.print();
};

//!It computes the DNA reverse
void convertN(char tmp[],const int& num);
//!It searchs in the hash table the kmerSizes of all reads, and return the pool contaning such window. 
//void Search(int l, int u);
//!It splits each read into the corresponding k-mers and store them in a hashtable. 
unsigned int readReads(ifstream& in, struct readFasta* reads);
int readOverlaps(int i, struct overlap* overlaps, ofstream& out);
void OverlapToDeconv(int id, int readNo, const SpMatCol& phi, int* bacPools, struct overlap* overlap, ofstream& out);

//void genErrors(int readNo, struct readFasta* read);
void FastaToKmers(int readNo, char* header, char* read, int pos, int p1, int p2);
//void showErrors(char* header, int* kmers, int numKmers, int p1, int p2);

void deconvOverlap(const SpMatCol& phi, int* bacPools, char* header, int numKmers, const SpMatCol& Y, ofstream& out);
void ReadToKmers(int id, int l, int u);
void OverlapToKmers(int id, int l, int u);
void OverlapToKmersLight(int id, int l, int u);
void RemoveLowComplexity(int id, int l, int u);
void ReadToKmersBack(int l, int u);
void CheckCorrect() {
	sh.checkCorrect();
}
void CheckCorrectBacs() {
	sh.checkCorrectBacs();
}

void PrintKmers(ostream& out) {
	sh.print(out);
}

void RemoveRepKmers() {
	cout << "[classIO:RemoveRep]Before removing rep kmers sh has " << sh.size() << " elements" << endl;
	sh.removeRepKmers();
	cout << "[classIO:RemoveRep]After removing rep kmers sh has " << sh.size() << " elements" << endl;
}

bool searchCopy(A_Type& b) {
	return sh.searchCopy(b);
}


bool searchUpdatePool(const A_Type& b, unsigned int pool, int val){
	return sh.searchUpdatePool(b, pool, val); 	
}

void computeProb() {
	sh.computeProb();
};

//!Performs read correction
//FIXME
//void correctReads(int pool, struct readFasta* reads, int numReads, HashTable<mysig>& pools2bacs);
//void dfs(int pool, ofstream& out, char*header, char* read, char* original, int direction, int pos, int kmerPos, unsigned short crtCset, int& prevECset, int& nextBCset, vector<struct cset>& csets, mybitsetx& e, mybitsetx& c, unsigned short orient, HashTable<mysig>& pools2bacs, unsigned int& numSearch, unsigned int& numHit, vector<int>& posChanged, bool& readCorrect);
//bool correct(int pool, char* read, int direction, int pos, int kmerPos, unsigned short crtCset, vector<struct cset>& csets, const mybitsetx& e, mybitsetx c, unsigned short orient, HashTable<mysig>& pools2bacs, unsigned int& numSearch, unsigned int& numHit);


//!It saves in a file the kmerSizes pools frequencies and their means.
void Average(vector<class OPERATOR>::iterator it,int i);
//!It saves in a file the kmerSizes corrections.
void Count(int i);

private:
	//!It encodes a buffer[num] and the pool id on two bitsets and  stores them into the hash table. It returns true if buffer cantains only [A,G,T,C].
	void  insert(const char buffer[],unsigned int num,unsigned int pool);
	void  insertLight(HashTable<mybitset>& sh, HashTable<mybitset> me[], HashTable<mybitset> others[],  			const char buffer[], int num, int pool);
	//!It encodes a buffer[num] in bitsets, and it searchs  into the hash table. It returns true if buffer cantains only [A,G,T,C].
	//inline mybitsetBAC search(class HashTable<mybitset>& rd,const char buffer[],const int& num,ofstream& out);
	//!It puts all the k-mers of the read stored in buffer in the hashtable rd. 
	//inline void readToKmersBack(class HashTable<mybitsetx>& rd,const char buffer[],const int& num);
	
//FIXME
/*friend void copyAll(const IOF<mybitsetx>& ref, IOF<mybitsetxbac>& ref1) 
{
	for (unsigned int i=0;i<POOLS;i++) {
		ref1.readCounts[i]=ref.readCounts[i];
		ref1.kmerCounts[i]=ref.kmerCounts[i];
	};
	//copyAllKmers(ref.sh, ref1.sh);
	//copyAllBacs(ref.shBAC, ref1.shBAC);
};*/

};
//! It uses for threads
struct thread_data{
   unsigned int l;
   unsigned int u;
};

//! It prints the read signature and the BAC list
//inline int printBitset(mybitsetBAC& I,ofstream& out);


class myhash
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


template <class A_Type>
void IOF<A_Type>::insertLight(HashTable<mybitset>& shLight, 
	HashTable<mybitset> me[F], HashTable<mybitset> others[F], const char buffer[], int num, int pool)
{
	
bool found;
int f;
class mybitset b,d,e;
for(unsigned int k=0; k<=num-(kmerSize)-1; k++)
	{
	int z=0;
	int j=DIM*2-1;
	bool findN=false;
	for (unsigned int c = 0; c<kmerSize; c++) 
		{
              	switch (buffer[c+k])
			{
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
	if (!findN)
	{
		kmerCounts[pool]++;
		f = F - 1;
		found = false;
		if (d<b) 
			e=d;
		else
			e=b;
		
		/*insert e into the hashtable*/
		
		//if (pool >0)
		//	cout << "[InsertLight]inserting " << e << " into the hashtable " << endl;

		if (shLight.search(e)) { 
				//if (pool >0)
				//	cout << "[insertLight]kmer found in sh" << endl;
				
				//if the k-mer is found in keep do nothing 
			found = true; 
		} 
		while ((found==false) && (f >= 0)) {
				
			if (others[f].search(e)) {
				
				//if (pool >0)
				//	cout << "[insertLight]kmer found in others[" << f << "]" << endl;
				
				found = true;
				others[f].removeCopy(e);
				
				if (f == (F - 1)) {//move the k-mer to keep 
						
					//if (pool>0)
					//	cout<< "[insertLight] pool " << pool<< " " << f << "==" << (F-1) << endl;
					
					shLight.insert(e);
				} 
				else { //move the k-mer from others[f] to me[f+1]
						
					//if (pool >0)
					//	cout << "[insertLight]pool  " << pool << f << "!=" << (F-1) << endl;
					
					me[f+1].insert(e);
				}
			} else if (me[f].search(e)) { 
				
				//if (pool >0)
				//	cout << "[insertLight]kmer found in me[" << f << "]" << endl;
				
				//if the k-mer is found in me do nothing 
				found = true; 
			};
			
			f--;

		} //end while ((!found) && (f >= 0))
			
		if (found == false) {
				
			//if (pool >0)
			//	cout << "[insertLight]kmer " << e << " about to be inserted in me[0]" << endl;
			
			me[0].insert(e);
		}; 	
	} //if (!findN)
	else
		numDiscard++;
	}
};


template <class A_Type>
void IOF<A_Type>::insert(const char buffer[],unsigned  int num,unsigned int pool){
A_Type b,d;
for(unsigned int k=0; k<=num-kmerSize-1; k++)
	{
	int z=0;
	int j=DIM*2-1;
	bool findN=false;
	for (unsigned int c = 0; c<kmerSize; c++) 
		{
              	switch (buffer[c+k])
			{
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
		kmerCounts[pool]++;
		if (d<b) 
			sh.searchIncPool(d,pool);
		else  
			sh.searchIncPool(b,pool);
	} else
		numDiscard++;
	}
};

template <class A_Type>
void IOF<A_Type>::ReadLight(){

using namespace general;

clock_t startGlobal,endGlobal;

/*Inits*/

HashTable<mybitset> shLight;
HashTable<mybitset> me[F];
HashTable<mybitset> others[F];
/*shLight.allocHash();
for (int f=0; f<F; f++) {
	me[f].allocHash();
	others[f].allocHash();
};*/
//FIXME
//shLight.allocHash(5);
shLight.allocHash(3);
me[0].allocHash(1);
me[1].allocHash(0);
//others[0].allocHash(12);
others[0].allocHash(8);
others[1].allocHash(2);

unsigned int num;
char buffer[MAXSIZE];

ifstream in;
//FILE* in;

char delimC[] = ".";
Parser parser;
parser.update(delimC,ReadFname);
/*string fileExt="";
if (parser.size()>1)
	{
	fileExt=parser.get(1);
	}
*/
//DEBUG
cout << endl;
cout << "[ReadLight]shLight size is " << shLight.size() << " and capacity " << shLight.capacity() << endl; 
for (int f =0; f<F; f++) {
	cout << "[ReadLight]me[" <<  f << "] size is " << me[f].size() << " and capacity " << me[f].capacity() << endl;
	cout << "[ReadLight]others[" <<  f << "] size is " << others[f].size() << " and capacity " << others[f].capacity() << endl;
	}
cout << endl;	


for (unsigned int i=0;i<numFiles;i++)
{
ostringstream of;
of<<parser.get(0)<<i<<"."<<parser.get(1);//<<"."<<parser.get(2)<<"."<<parser.get(3);

in.open(of.str().c_str(),ifstream::in);
//DEBUG
cout << endl << "***Pool " << of.str() << "***" << endl;
if(!in) 
	{
	cerr << "\nError opening input file " << of.str() << endl;
	exit(EXIT_FAILURE);
	}

startGlobal=clock();
while (in.good()) {
	buffer[0]='\0';
	//if (fgets(buffer, MAXSIZE, in)==NULL) 
	//	break;
	in.getline(buffer,MAXSIZE);
	//num=in.gcount();
	num=strlen(buffer);
	/*if (buffer[num-1]=='\n') {
		buffer[num-1]='\0';
		num--;
	}*/
	if ((num>=kmerSize+1)&&((buffer[0]=='A')||(buffer[0]=='G')||(buffer[0]=='C')||(buffer[0]=='T')||(buffer[0]=='N')))
		{
			readCounts[i]++;
			insertLight(shLight, me, others, buffer, num, i);
			
			/*if(readCounts[i]%500000==0) {
				cout<<"[ReadLight]Reads so far: "<<readCounts[i]<<", stored kmers: "<<shLight.size()<<", numDiscard "<<numDiscard<<endl;
			}*/	
		}
	};//end while(!in.eof())

numReads+=readCounts[i];

in.close();
//fclose(in);
endGlobal=clock();
cout << "\nTime to read pool "<<i<<": "<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl;
cout<<"Reads in this pool: "<<readCounts[i]<<" kmers: "<<kmerCounts[i]<<" kmers stored so far: "<<shLight.size() << " k-mers discarded " << numDiscard <<endl;

//DEBUG
cout << "\n[ReadLight]After processing pool " << i << endl;
cout << "shLight has size " << shLight.size() << ", capacity " << shLight.capacity() << " and ratio " << shLight.size()/((double)shLight.capacity()) << endl; 
for (int f = 0; f < F; f++) {
	cout << "me[" << f << "] has size " << me[f].size() << ", capacity " << me[f].capacity() << " and ratio " << me[f].size()/((double)me[f].capacity()) << endl; 
	cout << "others[" << f << "] has size " << others[f].size() << ", capacity " << others[f].capacity() << " and ratio " << others[f].size()/((double)others[f].capacity()) << endl; 
};
cout << endl;


//after processing all reads for pool i copy all me hashtables into the others hashtable 
//before continuing with the next pool
//if (i < (numFiles-1)) {
	for (int f = 0; f < F; f++) {
		//DEBUG
		cout << "[ReadLight]Moving me[" << f << "] into others[" << f << "]..." << endl;  
		others[f].insert(me[f]);
		cout << "[ReadLight]Done moving." << endl;
		me[f].makeEmpty();
	};
//};

//DEBUG
cout << "\n[ReadLight]After processing pool " << i << " and moving me to others " << endl;
cout << "shLight has size " << shLight.size() << ", capacity " << shLight.capacity() << " and ratio " << shLight.size()/((double)shLight.capacity()) << endl; 
for (int f = 0; f < F; f++) {
	cout << "me[" << f << "] has size " << me[f].size() << ", capacity " << me[f].capacity() << " and ratio " << me[f].size()/((double)me[f].capacity()) << endl; 
	cout << "others[" << f << "] has size " << others[f].size() << ", capacity " << others[f].capacity() << " and ratio " << others[f].size()/((double)others[f].capacity()) << endl; 
};
cout << endl;


}//end for each pool

/*cout << "[ReadLigh]Before calling makeEmpty() shLight has " << shLight.countElements() << " elemnts" << endl;
for (int f=0; f<F;f++) {
	cout << "[ReadLigh]Before calling makeEmpty() me[" << f << "] has " << me[f].countElements() << " elemnts" << endl;
	cout << "[ReadLigh]Before calling makeEmpty() others[" << f << "] has " << others[f].countElements() << " elemnts" << endl;
}*/

//delete all me and others hash tables
/*for (int f=0; f<F;f++) {
	me[f].makeEmpty();
	others[f].makeEmpty();
}*/

/*delete [] me;
delete [] others;*/

/*cout << "[ReadLigh]After calling makeEmpty() shLight has " << shLight.countElements() << " elemnts" << endl;
for (int f=0; f<F;f++) {
	cout << "[ReadLigh]After calling makeEmpty() me[" << f << "] has " << me[f].countElements() << " elements" << endl;
	cout << "[ReadLigh]After calling makeEmpty() others[" << f << "] has " << me[f].countElements() << " elemnts" << endl;
}*/


//resize down the hashtable if necessary
//FIXME
/*if (shLight.size()/((double)shLight.capacity()) <= MIN) {  
	cout <<"[ReadLight]shLight has size " << shLight.size() << ", capacity " << shLight.capacity() << " ratio " << shLight.size()/((double)shLight.capacity()) << " resizing down!!! " << endl;
	shLight.resize(0);
}*/

//copy shLight to sh
//cout << "[ReadLight]Inserting all the elements from shLight into sh..." << endl;
//insertAll(shLight, sh);
//cout << "[ReadLight]After insertion shLight has size " << shLight.size() << " and capacity " << shLight.capacity() << endl; 
//cout << "[ReadLight]After insertion sh has size " << sh.size() << " and capacity " << sh.capacity() << endl; 

//write shLight to disk instead to save space
shLight.Write(parser.get(0));

cout <<"\n[ReadLight]Done reading " << numFiles << " pools for a total of " << numReads << " reads." << endl;
};


template <class A_Type>
void IOF<A_Type>::ReadIntermediate() {
	using namespace general;
	
	char delimC[] = ".";
	Parser parser;
	parser.update(delimC,ReadFname);
	
	HashTable<mybitset> shLight;
	shLight.Read(parser.get(0));
	cout << "\n[ReadIntermediate]After reading from file shLight has size " << shLight.size() << " and capacity " << shLight.capacity() << endl; 
	
	cout << "[ReadIntermediate]Inserting all the elements from shLight into sh..." << endl;
	insertAll(shLight, sh);
	cout << "[ReadIntermediate]After insertion shLight has size " << shLight.size() << " and capacity " << shLight.capacity() << endl; 
	cout << "[ReadIntermediate]After insertion sh has size " << sh.size() << " and capacity " << sh.capacity() << endl; 

};


template <class A_Type>
void IOF<A_Type>::Read(){

using namespace general;

clock_t startGlobal,endGlobal;

unsigned int num;
char buffer[MAXSIZE];

ifstream in;
//FILE* in;

char delimC[] = ".";
Parser parser;
parser.update(delimC,ReadFname);

//DEBUG
cout << "[classIO:Read]Before processing all pools sh has size " << sh.size() << " and capacity " << sh.capacity() << endl; 

for (unsigned int i=0;i<numFiles;i++)
{
ostringstream of;
of<<parser.get(0)<<i<<"."<<parser.get(1);//<<"."<<parser.get(2)<<"."<<parser.get(3);
//DEBUG
cout << endl << "***Pool " << of.str() << "***" << endl;

in.open(of.str().c_str(),ifstream::in);
if(!in) 
	{
	cerr << "\nError opening input file " << of.str() << endl;
	exit(EXIT_FAILURE);
	}

startGlobal=clock();
while (in.good())
	{
	buffer[0]='\0';
	//if (fgets(buffer, MAXSIZE, in)==NULL)
	//	break;
	in.getline(buffer,MAXSIZE);
	//num=in.gcount();
	num=strlen(buffer);
	/*if (buffer[num-1]=='\n') {
		buffer[num-1]='\0';
		num--;
	}*/
	if ((num>=kmerSize+1)&&((buffer[0]=='A')||(buffer[0]=='G')||(buffer[0]=='C')||(buffer[0]=='T')||(buffer[0]=='N')))
		{
			insert(buffer,num,i);
			/*if(numReads%20000000==0) {
				cout<<"[Read]reads so far: "<<numReads<<" : "<<" stored kmers: "<<sh.size()<<" rejected  kmers: "<<numDiscard<<endl;
				//sh.print();
				sh.info();
			}	*/
		}
	}//end while (!in.eof())

in.close();
//fclose(in);
endGlobal=clock();
cout << "\n[classIO:Read]Time to read pool "<<i<<": "<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl;
//cout<<"Reads: "<<readCounts[i]<<" kmers: "<<kmerCounts[i]<<" kmers stored so far: "<<sh.size()<<endl;
};

//DEBUG
cout << "\n[classIO:Read]After reading all pools sh has size " << sh.size() << " and capacity " << sh.capacity() << endl << endl; 

};

template <class A_Type>
void IOF<A_Type>::ReadBac(){

using namespace general;

char delimC[] = "\t,";
int  PoolsBac[BACS][8];
Parser parser;
ifstream in;

in.open(BACFname.c_str(),ofstream::in);
if(!in) 
	{
	cerr << "\nError opening " << BACFname  << endl;
	exit(EXIT_FAILURE);
	}
char buffer[MAXSIZE];
int j=0;
while (!in.eof()) {
	buffer[0]='\0';
	in.getline(buffer,MAXSIZE);
	if (buffer[0]!='\0') {
		int num=in.gcount();
		if (buffer[num-1]!='\0') {
			buffer[num]='\0';
		}
		parser.update(delimC,buffer);	
		for (unsigned int i=0; i< parser.size(); i++ ) {
               		PoolsBac[j][i]=(atoi (parser.get(i).c_str()))-1;
		}
		j++;
	}
}	
in.close();
for (unsigned int j=0; j<BACS; j++ ) {
	mybitsetBAC p;
	p.set(PoolsBac[j][0]);//set BAC ID 
	for (unsigned int i=1; i<8; i++ ) {
               	p.setPool(PoolsBac[j][i]);
		//cout << PoolsBac[j][i] << " "; 
	}
	//cout << endl;
	if (shBAC.search(p)==0)
		shBAC.insert(p);
}

cout << "shBAC size " << shBAC.size() << endl;
};

/*void copyAll(const IOF<mybitsetx>& ref, IOF<mybitsetxbac>& ref1) 
{
	for (unsigned int i=0;i<POOLS;i++) {
		ref1.readCounts[i]=ref.readCounts[i];
		ref1.kmerCounts[i]=ref.kmerCounts[i];
	};
	//copyAllKmers(ref.sh, ref1.sh);
	//copyAllBacs(ref.shBAC, ref1.shBAC);
};*/


template <class A_Type>
unsigned int IOF<A_Type>::readReads(ifstream& in, struct readFasta* reads) {

unsigned int num;
char head[MAXSIZE];
char seq[MAXSIZE];

int index=-1;
struct readFasta* crt;
unsigned int numReads=0;

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
		//if (crt->header==NULL)
		//	crt->header=(char*)malloc(MAXSIZE*sizeof(char));
		//else
			crt->header[0]='\0';
		strncpy(crt->header, head, num);
		crt->header[num]='\0';
	
		//DEBUG
		//printf("%s\n", crt->header);
	
		seq[0]='\0';
		in.getline(seq, MAXSIZE);
		num=strlen(seq);
		
		assert((seq[0]=='A')||(seq[0]=='G')||(seq[0]=='C')||(seq[0]=='T')||(seq[0]=='N')||(seq[0]=='R'));
		//if (crt->read==NULL)
		//	crt->read=(char*)malloc(MAXSIZE*sizeof(char));
		//else
			crt->read[0]='\0';
		strncpy(crt->read, seq, num);
		crt->read[num]='\0';
			
	}//end if
}//end while

endGlobal=clock();
cout << "\n[read2kmers:readReads]Time to read pool: "<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl;
	
return numReads;

};


/*Builds hashtable on all reads (after correcting them using sga)
and resizes the hashtable if necessary. 
Avoids ReadLight() followed by Read() which is extremely time consuming*/
template <class A_Type>
void IOF<A_Type>::ReadBarley(){

using namespace general;
//using namespace std;
//using namespace MYBIT;

clock_t startGlobal,endGlobal;

unsigned int num;
char buffer[MAXSIZE];

ifstream in;
//FILE* in;

char delimC[] = ".";
Parser parser;
parser.update(delimC,ReadFname);

cout << "[ReadBarley]sh size is " << sh.size() << " and capacity " << sh.capacity() << endl; 
for (unsigned int i=0;i<numFiles;i++)
{
//DEBUG
cout << endl << "***Pool " << i << "***" << endl;

ostringstream of;
of<<parser.get(0)<<i<<"."<<parser.get(1);//<<"."<<parser.get(2)<<"."<<parser.get(3);

in.open(of.str().c_str(),ifstream::in);
//in=fopen(of.str().c_str(), "r");
/*if(in==NULL) 
	{
	cerr << "\nError opening input file " << of.str() << endl;
	exit(EXIT_FAILURE);
	}
*/
startGlobal=clock();
//while (!feof(in))
while (in.good())
	{
	buffer[0]='\0';
	//if (fgets(buffer, MAXSIZE, in)==NULL)
	//	break;
	in.getline(buffer,MAXSIZE);
	//num=in.gcount();
	num=strlen(buffer);
	/*if (buffer[num-1]=='\n') {
		buffer[num-1]='\0';
		num--;
	}*/
	if ((num>=kmerSize+1)&&((buffer[0]=='A')||(buffer[0]=='G')||(buffer[0]=='C')||(buffer[0]=='T')||(buffer[0]=='N')))
		{
			readCounts[i]++;	
			insert(buffer,num,i);
			/*if(numReads%20000000==0) {
				cout<<"[Read]reads so far: "<<numReads<<" : "<<" stored kmers: "<<sh.size()<<" rejected  kmers: "<<numDiscard<<endl;
				//sh.print();
				sh.info();
			}	*/
		}
	}//end while (!in.eof())

numReads+=readCounts[i];

//fclose(in);
in.close();

endGlobal=clock();
cout << "\n[classIO:ReadBarley]Time to read pool "<<i<<": "<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl;
cout<<"[classIO:ReadBarley]Reads in this pool: "<<readCounts[i]<<", kmers: "<<kmerCounts[i]<<" kmers stored "<<sh.size() << " k-mers discarded " << numDiscard <<endl;

}//end for i

//DEBUG
cout <<"\n[classIO:ReadBarley]Done reading " << numFiles << " pools for a total of " << numReads << " reads." << endl;
cout << "\n[classIO:ReadBarley]After reading all pools sh has size " << sh.size() << " and capacity " << sh.capacity() << endl << endl; 
sh.checkCorrect();

};

};
