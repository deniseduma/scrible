/***************************************************************************
 *   Copyright (C) 2011 by Marco Beccuti   *
 *   beccuti@di.unito.it   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef __CON_H__
	#define __CON_H__
	#include "const.h"
#endif

#ifndef __FSTREAM__
	#define __FSTREAM__
	#include <fstream>
#endif

#ifndef __IOS_H__
	#define __IOS_H__
	#include <iostream>
#endif 


#ifndef __STDL__
	#define __STDL__
	#include <stdlib.h>
#endif

#ifndef __VCT_H__
	#define __VCT_H__
	#include <vector>
#endif

#ifndef __CMB_H__
	#define __CMB_H__
	#include "classMybitset.h"
#endif


#ifndef __CMBX_H__
	#define __CMBX_H__
	#include "classMybitsetX.h"
#endif

#ifndef __CMBXB_H__
	#define __CMBXB_H__
	#include "classMybitsetXBac.h"
#endif

#ifndef __CMS_H__
	#define __CMS_H__
	#include "classMySig.h"
#endif

#ifndef __CMBO_H__
	#define __CMBO_H__
	#include "classMybitsetChar.h"
#endif

#ifndef __LMT_H__
	#define __LMT_H__
	#include <limits.h>
#endif


#ifndef __OPER_H__
	#define __OPER_H__
	#include "operation.h"
#endif

#include "trie.h"

#include <omp.h>
#include <assert.h>
#include <typeinfo>
#include <sstream>

#include <typeinfo>
#include <math.h>

extern unsigned int hashSizes[SIZES];

namespace Cl_HASH
{
using namespace  MYBIT;
using namespace std;
using namespace OPER;



template <class A_Type> 
class HashNode {
public: 
    A_Type b;
    HashNode *next;
    HashNode(A_Type& bb, HashNode<A_Type>* n): b(bb), next(n) {
	//b=bb;
	//next=n;
    }
}; 


template <class A_Type> 
class HashTable {
private:
	typedef HashNode<A_Type>* PNode;

	PNode* HashVec;
	//!Hash table size, an index into hashSizes (see above extern array)  
	unsigned short SizeHash;
	unsigned int cpty;
	//!Number of elements in the hash table
	unsigned int NumEl;

public:

//!Empty Constructor.
HashTable(){
	NumEl=0;
	SizeHash=0;
	cpty=0;
	HashVec=NULL;
	//HashVec=new PNode[hashSizes[SizeHash]];
	//for (unsigned int i=0;i<hashSizes[SizeHash];i++)
	//	HashVec[i]=NULL;
	//cout << "Constructor of hashtable called." << endl;
};
	
	
//!Destructor
~HashTable(){
	makeEmpty();
	delete [] HashVec;
	HashVec=NULL;
	//DEBUG
	cout << "Destructor of hashtable called." << endl;
};

//!It updates the hash size
void allocHash() {
	//SizeHash++;
	cpty=hashSizes[SizeHash];
	HashVec=new HashNode<A_Type>*[cpty];
	//HashVec = (PNode*)malloc(SizeHash*sizeof(PNode));
	for (unsigned int i=0;i<cpty;i++)
		HashVec[i]=NULL;
};

//!It updates the hash size
void allocHash(unsigned short sizeHash) {
	//SizeHash++;
	SizeHash=sizeHash;
	cpty=hashSizes[SizeHash];
	HashVec=new HashNode<A_Type>*[cpty];
	//HashVec = (PNode*)malloc(SizeHash*sizeof(PNode));
	for (unsigned int i=0;i<cpty;i++)
		HashVec[i]=NULL;
};

//!It resizes the hashtable 
bool resize(short up){
	
	/*save old hash table*/
	unsigned int oldEl=NumEl;
	unsigned short oldSize = SizeHash;
	unsigned int oldCpty=hashSizes[oldSize];
	PNode* oldHashVec = HashVec;
	
	if (up==1) { 	
		if (SizeHash < SIZES) {
			SizeHash++;
			cpty=hashSizes[SizeHash];
			HashVec=new HashNode<A_Type>*[cpty];
		}	
		else
			return false;
	}		
	else if (up==0) { 	
		if (SizeHash > 0) {
			SizeHash--;
			cpty=hashSizes[SizeHash];
			HashVec=new HashNode<A_Type>*[cpty];
		}	
		else
			return false; 
	}		
	
	NumEl=0;
	for (unsigned int i=0; i<cpty; i++)
		HashVec[i]=NULL;
	
	/*rehash all the elements from the old hashtable into the new one*/
	for(unsigned int i=0; i<oldCpty; i++) {
		PNode tmpHashNode = oldHashVec[i];
    		while(tmpHashNode != NULL) {
			insert(tmpHashNode->b);
			tmpHashNode = tmpHashNode->next;
    		}
	}
	cout << "[resize]oldCapacity " << oldCpty << " new capacity " << cpty << endl;
	cout << "[resize]oldEl " << oldEl << " NumEl " << NumEl << " old ratio " << (oldEl/((double)oldCpty)) << " new ratio " << (NumEl/(double)cpty) << endl;
	cout << "[resize]oldSize " << oldSize << ", SizeHash " << SizeHash << endl;
	
	assert(oldEl==NumEl);
	
	/*empty the old hashtable*/

	HashNode<A_Type>* tmpHashNode;
	for (unsigned int i =0; i<oldCpty; i++) {
		while (oldHashVec[i] != NULL) {
			tmpHashNode = oldHashVec[i];
			oldHashVec[i] = oldHashVec[i]->next;
			delete tmpHashNode;
			oldEl--;
		}
		
	}
	assert(oldEl==0);
	
	delete [] oldHashVec;

	return true;
};


//!It inserts the input value in the hash table, using chaining.  
inline void insert(A_Type& b) {
		
	unsigned int probe = hash(b.trans());
	assert(probe < capacity());
	HashVec[probe] = new HashNode<A_Type>(b, HashVec[probe]);
	NumEl++;
	
	if (size() / ((double)capacity()) > MAX) {
		//DEBUG
		cout << "[hashChaining:insert]Hash table needs resizing." << endl;
		if (!resize(1)) { 
			cout << "Error resizing hash table up." << endl;
			exit(1);
		}	
		//DEBUG
		cout << "[hashChaining:insert]Resized crt hash table up to " << capacity() << endl;
	}
};


//!Counts the elemnts in the hashtable
inline unsigned int countElements() {
	unsigned int count =0;
	for (unsigned int i=0; i<capacity(); i++) {
		PNode tmpHashNode=HashVec[i];
		while (tmpHashNode != NULL) {
			count++;
			tmpHashNode=tmpHashNode->next;
		}
	}
	assert(count==NumEl);

	return count;
}

//!Hash function. It takes in input a integer value.
inline unsigned int hash(const unsigned long long& value) const{
	//return (H(value)) % SizeHash;
	register unsigned  long long u=value;
	register unsigned  long long v=value>>32;
	//return (((u*2+v)*A)%L)%SizeHash;
	return ((u*2+v)*A)%capacity();
	//return (value%L) % SizeHash;
};

//!It inserts all the elements from the input hash table into the current hash table. 
void insert(const HashTable<A_Type>& sh){
	for (unsigned int i=0;i<sh.capacity();i++)
		{
		PNode tmpHashNode = sh.HashVec[i];
    		while(tmpHashNode != NULL)
		{
		 	insert(tmpHashNode->b);
        		tmpHashNode = tmpHashNode->next;
    		}
	}
};		

//!Checks that all kmers stored in the hashtable appear in at least F pools
void  checkCorrect() {
	unsigned int freqs[POOLS]; 
	for (int i=0; i<POOLS; i++)
		freqs[i]=0;
	for (unsigned int i=0;i<capacity();i++) {
		PNode tmpHashNode = HashVec[i];
		while (tmpHashNode != NULL) {
			//assert(tmpHashNode->b.count()>=MINREP);// && tmpHashNode->b.count()<=MAXREP);
			freqs[tmpHashNode->b.getNumPools()-1]++;
			tmpHashNode = tmpHashNode->next;
		}
	}
	cout << "[hashChaining:Check correct]" << endl;
	for (int i=0;i<POOLS;i++)
		cout<<(i+1)<<" "<<freqs[i]<<endl; 
};


//!Checks that all kmers stored in the hashtable appear in at least F pools
void  checkCorrectBacs() {
	unsigned int freqs[POOLS]; 
	unsigned int freqsBacs[BACS];
	for (unsigned int i=0; i<POOLS; i++)
		freqs[i]=0;
	for (unsigned int i=0; i<BACS; i++)
		freqsBacs[i]=0;
	for (unsigned int i=0;i<capacity();i++) {
		PNode tmpHashNode = HashVec[i];
		while (tmpHashNode != NULL) {
			//assert(tmpHashNode->b.count()>=MINREP);// && tmpHashNode->b.count()<=MAXREP);
			unsigned short numPools=tmpHashNode->b.getNumPools();
			freqs[numPools-1]++;
			if (numPools>=LOW2&&numPools<=HIGH2)  
				freqsBacs[tmpHashNode->b.getNumBacs()-1]++;
			tmpHashNode = tmpHashNode->next;
		}
	}
	cout << "[hashChaining:Check correct bacs]" << endl;
	for (unsigned int i=0;i<POOLS;i++)
		cout<<(i+1)<<" "<<freqs[i]<<endl; 
	for (unsigned int i=0;i<BACS;i++)
		cout<<(i+1)<<" "<<freqsBacs[i]<<endl; 
};


//!It inserts all the elements from the first hash table into the second hash table. 
//template<class A_Type2>
friend void insertAll(const HashTable<mybitset>& sh1, HashTable<A_Type>& sh2) {
	for (unsigned int i=0;i<sh1.capacity();i++)
	{
		HashNode<mybitset>* tmpHashNode = sh1.HashVec[i];
		while(tmpHashNode != NULL) {
		 	//mybitsetx mb(tmpHashNode->b.getBitset());
		 	A_Type mb(tmpHashNode->b);
			sh2.insert(mb);
			tmpHashNode = tmpHashNode->next;
    		};
	};
	assert(sh2.size()==sh1.size());
};


/*friend void copyAllKmers(const HashTable<mybitsetx>& sh, HashTable<mybitsetxbac>& sh1) {
	for (unsigned int i=0;i<sh.capacity();i++)
	{
		HashNode<mybitsetx>* tmpHashNode = sh.HashVec[i];
		while(tmpHashNode != NULL) {
		 	//mybitsetx mb(tmpHashNode->b.getBitset());
		 	mybitsetxbac mb(tmpHashNode->b);
			sh1.insert(mb);
			tmpHashNode = tmpHashNode->next;
    		};
	};
	assert(sh1.size()==sh.size());
};


friend void copyAllBacs(const HashTable<mybitsetBAC>& sh, HashTable<mybitsetBAC>& sh1) {
	for (unsigned int i=0;i<sh.capacity();i++)
	{
		HashNode<mybitsetBAC>* tmpHashNode = sh.HashVec[i];
		while(tmpHashNode != NULL) {
		 	//mybitsetx mb(tmpHashNode->b.getBitset());
		 	mybitsetBAC mb(tmpHashNode->b);
			sh1.insert(mb);
			tmpHashNode = tmpHashNode->next;
    		};
	};
	assert(sh1.size()==sh.size());
};*/


//template<class A_Type2>
/*void cmpTables(const HashTable<mybitset>& sh) {
	for  (unsigned int i=0;i<sh.capacity();i++) 
	{
		HashNode<mybitset>* tmpHashNode = sh.HashVec[i];
		while (tmpHashNode != NULL) {
			A_Type mb(tmpHashNode->b);
			assert(search(mb)==false);
			tmpHashNode = tmpHashNode->next;
		}
	}
};*/

//It empties the hash table
inline void makeEmpty() {
	//DEBUG
	//cout << "makeEmpty called" << endl;	
	HashNode<A_Type>* tmpHashNode;
	for (unsigned int i =0; i<capacity(); i++) {
		while (HashVec[i] != NULL) {
			tmpHashNode = HashVec[i];
			HashVec[i] = HashVec[i]->next;
			delete tmpHashNode;
			NumEl--;
		}
	}
	//DEBUG
	assert(NumEl==0);
	for (unsigned int i =0; i<capacity(); i++) 
		assert(HashVec[i]==NULL);
}


//!Search a value in the hash table  and if found copy in b the pool bitvector
inline bool searchCopy(A_Type& b) const{
	//PNode tmpHashNode = HashVec[hash(b.trans())];
	unsigned int probe = hash(b.trans());
	assert(probe < capacity());
	PNode tmpHashNode = HashVec[probe];
	while(tmpHashNode != NULL) {
		if(tmpHashNode->b == b) { 
			#pragma omp critical
			tmpHashNode->b.copy(b);
			return true;
		}	
		tmpHashNode = tmpHashNode->next;
    	}
	return false;
};
	
	
//!Search  a value in the hash table and if found return true.
inline bool search(A_Type& b) const {
	PNode tmpHashNode = HashVec[hash(b.trans())];
	while(tmpHashNode != NULL) {
       		if(tmpHashNode->b == b) {
			return true;
		}
       		tmpHashNode = tmpHashNode->next;
	}
return false;
};

//FIXME
//Has to be changed when running for rice or BARLEY. See commented code! 
//!Search  a value in the hash table and if found increment pool count.
inline bool searchIncPool(A_Type& b, int pool){
	
	PNode tmpHashNode = HashVec[hash(b.trans())];
	while(tmpHashNode != NULL) {
        	if(tmpHashNode->b == b) {
			tmpHashNode->b.incPool(pool);
			return true;
		}
        	tmpHashNode = tmpHashNode->next;
    	};
	//for barley
	//b.incPool(pool);
	//insert(b);
	//return true;
	return false;
};

//!Search  a value in the hash table and if found update pool count by val.
inline bool searchUpdatePool(const A_Type& b, unsigned int pool, short val){
	
	//bool found = false;
	PNode tmpHashNode = HashVec[hash(b.trans())];
	while(tmpHashNode != NULL) {
        	if(tmpHashNode->b == b) {
			#pragma omp critical 
			{
			tmpHashNode->b.updatePool(pool, val);
			}
			return true;
			//found = true;
			//break;
		}
        	tmpHashNode = tmpHashNode->next;
    	}
	/*if (!found) {
		b.updatePool(pool, val);
		insert(b);
	}*/
	return false;
};


//!Removes b from the hashtable.
bool removeCopy(A_Type& b){
	
	bool res = false;
	
	//cout << "[removeCopy]removeCopy called" << endl;
	unsigned int probe = hash(b.trans());
	if (HashVec[probe] == NULL) {
		return false;
	} 
	PNode tmpHashNode; 
	if (HashVec[probe]->b == b) {
		tmpHashNode = HashVec[probe];
		//tmpHashNode->b.copy(b);
		HashVec[probe] = HashVec[probe]->next;
		delete tmpHashNode;
		NumEl--;
		res=true; 	
		//DEBUG
		//cout << "[removeCopy] node found in list head"<< endl;
	} else {
	
		PNode prevHashNode = HashVec[probe];
		tmpHashNode=HashVec[probe]->next; 
		while(tmpHashNode != NULL) {
			if (tmpHashNode->b == b) 
				break;
       			prevHashNode = tmpHashNode;
			tmpHashNode=tmpHashNode->next;
		}	
		
		if (tmpHashNode != NULL) {
			//tmpHashNode->b.copy(b);
			prevHashNode->next = tmpHashNode->next;
			delete tmpHashNode;
			NumEl--;
			res=true;
		}
	}	
		
	/*resize hash table down if necessary*/
		
	/*if (res && SizeHash>0) {  
		if (size() / ((double)capacity()) < MIN) {
			resize(0);	
			cout << "[hashChaining:removeCopy]Hash table resized down to " << capacity() << " new SizeHash is " << SizeHash << " new ratio is " << size()/ ((double)capacity()) << endl;
		}
	}*/
		
	return res;
};


//!Search  a value in the hash table if found update the pool.
void removeRepKmers(){
	for (unsigned i=0; i<capacity();i++) {
		PNode tmpHashNode = HashVec[i];
		while(tmpHashNode != NULL) {
        		//if(tmpHashNode->b.count()<MINREP || tmpHashNode->b.count()>MAXREP) 
        		if(tmpHashNode->b.getNumPools()<MINREP) 
				removeCopy(tmpHashNode->b);
			 tmpHashNode = tmpHashNode->next;
    		};
	};
};


//!Print the linked list stored in HashTable[row] in a file
void printList(int row, ostream& out){
	PNode tmpHashNode = HashVec[row];
    	while(tmpHashNode != NULL) {
		out << tmpHashNode->b;
		tmpHashNode = tmpHashNode->next;
    	}
};


//!Print the linked list stored in HashTable[row]
void printList(int row){
	PNode tmpHashNode = HashVec[row];
	if (tmpHashNode != NULL)
		cout<<"\nHASH["<<row<<"]:\t"; 
    	while(tmpHashNode != NULL) {
		cout << tmpHashNode->b << " ";
        	tmpHashNode = tmpHashNode->next;
    	}
};

//!Print the size of the linked list stored in HashTable[row]
void printSize(int row){
	PNode tmpHashNode = HashVec[row];
	int i=0;
	if (tmpHashNode != NULL)
		cout<<"\nHASH["<<row<<"]:\t"; 
    	while(tmpHashNode != NULL) {
		i++;
        	tmpHashNode = tmpHashNode->next;
    	}
	if (i>0)
		cout<<i<<endl;
};

//!Print the whole HashTable in a file 
void print(ostream& out){
	for(unsigned int i=0;i<capacity();i++) {
        	printList(i,out);
	}
};


//!Print the whole HashTable 
void print(){
	for(unsigned int i=0;i<capacity();i++) {
        	printList(i);
	}
	cout << endl;
};

//!Return the HASH size
unsigned int size() const {
	return NumEl;
};

//!Return the capacity of the hashtable
unsigned int capacity() const {
	return cpty;
};


//!Print some statistical information on the Hash Table
void info() {
	unsigned int hashentry=0;
	for (unsigned int i=0;i<capacity();i++) {
		if (HashVec[i] != NULL)
			hashentry++;	
	}
	cout<<"Hash table buckets: "<<hashentry<<endl;
	cout<<"Avergage size of collisions list: "<<(double)NumEl/hashentry<<endl;
};
	
//!Take in input the output file name and save the hash table in binary format.
void Write(string output) {	
	
	FILE * pFile;
	PNode tmpHashNode;
	string tmp=output+".hash";
	pFile = fopen ( tmp.c_str(), "w" );
	if (pFile==NULL) {
		cerr<<"\nError opening output file " << tmp << endl;
		exit(EXIT_FAILURE);
	}
	
	//write hash size
  	fwrite (&SizeHash, sizeof(SizeHash),1 , pFile );
	//write number of elements
	fwrite (&NumEl, sizeof(NumEl), 1, pFile );
	//write the hash table elements
	for (unsigned int i = 0; i<capacity(); i++) {
		tmpHashNode = HashVec[i];
		while (tmpHashNode != NULL) {
			fwrite (&tmpHashNode->b, sizeof(A_Type), 1, pFile);
			//DEBUG
			//cout << "[hashChaining:Write]written b is " << tmpHashNode->b << endl;
			tmpHashNode = tmpHashNode->next;
		}
	}
	fclose (pFile);
	cout<<"[hashChaining:Write]Hash table saved in: "<<tmp<<endl;
};



//!Takes in input the intput file name and updates the hash table
bool Read(string output) {
	
	FILE * pFile;
	string tmp=output+".hash";
  	pFile = fopen ( tmp.c_str(), "r" );
	if (pFile == NULL) {
		cerr<<"\n[hashChaining:Read]Error opening input file " << tmp << endl;
		cout<<"\nHash Table binary encoding not found!"<<endl;
		return false;
	}
	
	//read hash size
  	if (fread (&SizeHash, sizeof(SizeHash),1 , pFile )!=1) {	
		cerr<<"Error reading hash size"<<endl; 
		exit (EXIT_FAILURE);
	}
	cout << "\n[hashChaining:Read]hash size read from file " << SizeHash << endl;
	
	unsigned int NumElements;
	//read number of elements
	if (fread (&NumElements, sizeof(NumEl), 1, pFile )!=1) {
		cerr<<"Error reading the number of elements"<<endl; 
		exit (EXIT_FAILURE);
	}
	cout << "[hashChaining:Read]NumElements read from file " << NumElements << endl;
	
	//FIXME
	allocHash(SizeHash);

	//read elements and insert them one by one into the hashtable
	for (unsigned int i = 0; i<NumElements; i++) {
		A_Type b;
		if (fread(&b, sizeof(A_Type), 1, pFile) != 1) {
			cerr << "Error reading element " << i<< endl; 
			exit (EXIT_FAILURE);
		} 
		insert(b);
	}

	assert(NumEl==NumElements);

	fclose (pFile);
	cout<<"[hashChaining:Read]Hash Table binary encoding found. Size is " << size() << " and capacity " << capacity() << endl;

	return true;
};

void RetrieveBacs(Trie* trie) {
	unsigned short numBacs, numPools; 
	unsigned short bacs[MAX_BACS];
	unsigned short pools[POOLS];
	for(unsigned int i=0;i<capacity();i++) {
		PNode tmpHashNode = HashVec[i];
    		while(tmpHashNode != NULL) {
			numPools=tmpHashNode->b.getPools(pools);
			assert(numPools>=3&&numPools<=91);
			if (numPools>=LOW2&&numPools<=HIGH2) {
				numBacs=0;
				trie->search(pools, numPools, bacs, numBacs);
				if (!(numBacs>=0&&numBacs<=2197)) {
					cout << "[RetrieveBacs]numBacs is " << numBacs << "!!!" << endl;
					exit(1);
				}
				tmpHashNode->b.setBacs(bacs, numBacs);
			}
			tmpHashNode = tmpHashNode->next;
    		}
	}
};

void computeProb() {
	PNode tmpHashNode;
	unsigned short numPools;
	for(unsigned int i=0;i<capacity();i++) {
		tmpHashNode = HashVec[i];
    		while(tmpHashNode != NULL) {
			numPools=tmpHashNode->b.getNumPools();
			
			//V1
			double avg1=0;
			for (unsigned int i=0;i<POOLS;i++) 
				avg1 += tmpHashNode->b.getPoolCount(i);
			
			avg1 /= numPools;
			avg1 -= 3;
			double prob1=1.0 / (1.0 + exp((-1)*avg1));
			
			//V2
			double avg2=1;
			for (unsigned int i=0;i<POOLS;i++) 
				if (tmpHashNode->b.getPoolCount(i)>0)
					avg2 *= tmpHashNode->b.getPoolCount(i);
			
			double prob2=1.0/pow(avg2, 1.0/numPools);

			cout << tmpHashNode->b << endl;
			cout << "probability 1 " << prob1 << endl;
			cout << "probability 2 " << prob2 << endl;
			
			tmpHashNode = tmpHashNode->next;
		};	
	};
};	


//!It takes in input the intput file name and updates the hash table
bool ReadOld(string output) {
	
	FILE * pFile;
	string tmp=output+".hash.old";
  	pFile = fopen ( tmp.c_str(), "r" );
	if (pFile==NULL)
  		{
		cerr<<"\n[hashChaining:ReadOld]Error opening input file " << tmp << endl;
		cout<<"\nHash Table binary encoding not found!"<<endl;
		return false;
		}
	cout << "hash file is " << tmp << endl;
	//write hash size
  	unsigned int SizeVector;
	if (fread (&SizeVector, sizeof(SizeVector),1 , pFile )!=1)	
		{	
		cerr<<"[ReadOld]Reading error"<<endl; 
		exit (EXIT_FAILURE);
		}
	cout << "[hashChaining:ReadOld]SizeVector is " << SizeVector << endl;
	unsigned int * HashVector=(unsigned int*)malloc(SizeVector*sizeof(unsigned int));	
	//for (unsigned int i=0;i<SizeVector;i++)
	//	HashVector[i]=DEFAULTP;
	//read int vector
	if (fread (HashVector, sizeof(unsigned int), SizeVector, pFile )!=SizeVector)
		{	
		cerr<<"[ReadOld]Reading error"<<endl; 
		exit (EXIT_FAILURE);
		}	
  	fclose (pFile);
	
	//DEBUG
	cout << "[hashChaining:ReadOld]HashVector succesfully read from file." << endl;
	
	tmp=output+".collist.old";
	pFile = fopen ( tmp.c_str(), "r" );
	if (pFile==NULL)
  		{
		cerr<<"\n*****[hashChaining:ReadOld]Error opening input file " << tmp << endl;
		cout<<"\nHash Table binary encoding not found!"<<endl;
		return false;
		}
	cout << "collist file is " << tmp << endl;
	//read collisions list size
	unsigned int SizeList;
	if (fread (&SizeList, sizeof(SizeList), 1, pFile )!=1)
		{
		cerr<<"[ReadOld]Error reading SizeList "<<endl; 
		exit (EXIT_FAILURE);
		}
	cout << "[hashChaining:ReadOld]SizeList is " << SizeList << endl;
	//read number of elements
	unsigned int NumElements;
	if (fread (&NumElements, sizeof(NumElements), 1, pFile )!=1)
		{
		cerr<<"[ReadOld]Error reading NumElements "<<endl; 
		exit (EXIT_FAILURE);
		}	
	cout << "[hashChaining:ReadOld]NumElements is " << NumElements << endl;
	//read myclassbit vector
	//mybitsetChar* ListCol=(mybitsetChar*)malloc(SizeList*sizeof(mybitsetChar));
	A_Type* ListCol=(A_Type*)malloc(SizeList*sizeof(A_Type));
	if (fread (ListCol, sizeof(A_Type), SizeList, pFile )!=SizeList)
		{
		cerr<<"[ReadOld]Error reading ListCol "<<endl; 
		exit (EXIT_FAILURE);
		}
	fclose (pFile);
	
	//DEBUG
	cout << "[hashChaining:ReadOld]ListCol succesfully read from file." << endl;
	
	/*SizeHash=0;
	HashVec=new HashNode<A_Type>*[hashSizes[SizeHash]];
	for (unsigned int i=0;i<hashSizes[SizeHash];i++)
		HashVec[i]=NULL;
	*/
	/*for (unsigned int i=0;i<SizeVector;i++) {
		if (HashVector[i]>SizeList && HashVector[i]!=DEFAULTP) 
			cout << "i=" << i << " " << HashVector[i] << endl;
		//assert(HashVector[i]==DEFAULTP || (HashVector[i]>=0 && HashVector[i]<SizeList));
	}*/
	/*for (unsigned int i=0; i<SizeList; i++) {
		if (ListCol[i].Inext>SizeList && ListCol[i].Inext!=DEFAULTP) 
			cout << ListCol[i] << endl << ListCol[i].Inext << endl;
		//cout << ListCol[i].getBitset() << " " << ListCol[i].getNumPools() << endl;
	}*/
	unsigned int tmpInode;
	for (unsigned int i=0; i<SizeVector; i++) {
		tmpInode=HashVector[i];
		while (tmpInode != DEFAULTP) {
			//if (ListCol[tmpInode].count() <= F) {
				//DEBUG
				//cout << "[hashChaining:ReadOld]" << b_old.getBitset() << " appears in " << b_old.getNumPools() << "pools" << endl;
				//insert(mybitset(ListCol[tmpInode].getBitset()));
			//}
			tmpInode=ListCol[tmpInode].Inext;
			if (tmpInode>SizeList && tmpInode != DEFAULTP){
				cout << "[" << i << "] " << tmpInode << endl;
				//exit(1);
			}
		}
	}
	
	//DEBUG
	cout << "[hashChaining:ReadOld]Old hash table succesfully created." << endl;
	
	//assert(NumEl==NumElements);

	free(ListCol);
	free(HashVector);

	return true;
};


};

};


