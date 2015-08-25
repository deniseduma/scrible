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


#ifndef __BIT_H__
	#define __BIT_H__
	#include <bitset>

#endif
#ifndef __MTH_H__
	#define __MTH_H__
	#include <math.h>
#endif

#ifndef __SET_H__
	#define __SET_H__
	#include <list>
#endif

#ifndef __CON_H__
	#define __CON_H__
	#include "const.h"
#endif

#include <cstring>

extern unsigned char c[8];
extern unsigned char r[8];
 
namespace MYBIT{

using namespace std;

class mybitset{

protected:
	//!Sequence  in binary format.
	bitset <DIM*2> cc;
public:


//!Empty Constructor.
mybitset(): cc() {
};

//!Constructor
mybitset(const bitset<2*DIM>& b): cc(b) {}

//!Assignmenet operator
mybitset& operator=(const mybitset& b) {
	if (this!=&b) {
		cc=b.cc;
	}	
	return *this;
}

//!Copy constructor
mybitset(const mybitset& b): cc(b.cc) {}

virtual ~mybitset() {}
	
//!It sets to val  the ith position  of the sequence's bitvector. 
inline void set(const int& i,int val){
	cc[i]=val;
};

//!It returns the value of ith position of the sequence's bitvector.
inline int get(const int& i){
	return cc[i];
};

//!It flips the bit at position i
inline void flip(const int& i) {
	cc.flip(i);
}

//!It flips the bit at position i
inline char getBase(const int& i) const{
	if (cc[2*i]==0)
		if (cc[2*i+1]==0)
			return 'A';
		else 
			return 'C';	
	else 
		if (cc[2*i+1]==0)
			return 'G';
		else 
			return 'T';
		
}

//!It flips the bit at position i
inline void setBase(const int& i, const char& base) {
	switch (base) {
	case 'A': 
		cc[2*i]=0;
		cc[2*i+1]=0;
		break;
	case 'C': 
		cc[2*i]=0;
		cc[2*i+1]=1;
		break;
	case 'G': 
		cc[2*i]=1;
		cc[2*i+1]=0;
		break;
	case 'T': 
		cc[2*i]=1;
		cc[2*i+1]=1;
		break;
	}	
}

inline static char rev(const char& base) {
	switch (base) {
		case 'A':
			return 'T';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
		case 'T':
			return 'A';
		default: 
			return 'N';
	}
}

//!It returns the bitset
inline bitset<2*DIM> getBitset() {
	return cc;
};

virtual inline void copy(mybitset& p) {};

//!It inverts a window encoded in the sequence's bitvector.
inline void invert(){
	bitset <DIM*2> inv;
	int i=0;
	for (int jj=DIM*2-1;jj>=0;jj=jj-2)
		{
   		if (cc[jj]==cc[jj-1])
                	if(cc[jj]==1)
                        	{
                        	inv[i]=inv[i+1]=0;
                        	}
                      	else
                        	{
                        	inv[i]=inv[i+1]=1;
                        	}
                else
                    	{
                    	inv[i]=cc[jj];
                   	 inv[i+1]=cc[jj-1];
                    	}
                   i=i+2;

		}
	cc=inv;
};


//!It compares the two sequences  p1 and p2  and returns 0 when p1=p2, 1 when p1>p2 and 2 when p1<p2. 
inline int compare(const mybitset& p2){
	int i=0;
 	while (((this->cc[i]==p2.cc[i]))&&(i<=DIM*2))
 			{
 			i++;		
 			}
 	if (i>DIM*2)
 			{
 			return 0;
 			}
 	if (this->cc[i]==1)
 			return 1;
 		else		
 			return 2;
};
	
//!Operator==
friend bool operator==(const mybitset& p1, const mybitset& p2){
	if (p1.cc==p2.cc) {
		return true;
	};
	return false;
};


//!Operato!=
friend bool operator!=(mybitset& p1,mybitset& p2){
	if (p1.cc!=p2.cc) {
		return true;
	};
	return false;
};


//!Operator< 
friend bool operator<(const mybitset& p1, const mybitset& p2){
  	int i=0;
  	while (((p1.cc[i]==p2.cc[i]))&&(i<=DIM*2))
		{
		i++;
		}
	if (i>DIM*2)
		{
		return false;
		}
	if (p1.cc[i]==1)
		return false;
	else		
		return true;
};


//!Operator<<
friend ostream& operator<<(ostream& out, const mybitset& p){
	p.print(out);
	return out;
}

virtual void print(ostream& out) const { 
	int zz=0;
	string buffer="";
	zz=0;
	for (int ii=0;ii<DIM*2;ii++)
		{
		if ((cc[ii]==0)&& (cc[ii+1]==0))
				buffer+='A';
		else
			if ((cc[ii]==0)&&(cc[ii+1]==1))
				buffer+='C';
			else
				if ((cc[ii]==1)&&(cc[ii+1]==0))
					buffer+='G';
				else
					buffer+='T';
			ii++;
			zz++;
		}
	
	out<<buffer<<endl;
};

string seq() const { 
	int zz=0;
	string buffer="";
	zz=0;
	for (int ii=0;ii<DIM*2;ii++)
		{
		if ((cc[ii]==0)&& (cc[ii+1]==0))
				buffer+='A';
		else
			if ((cc[ii]==0)&&(cc[ii+1]==1))
				buffer+='C';
			else
				if ((cc[ii]==1)&&(cc[ii+1]==0))
					buffer+='G';
				else
					buffer+='T';
			ii++;
			zz++;
		}
	
	return buffer;
};

//!It returns the sequence (buffer). 
inline void bit2buffer(char buffer[]){
	int zz=0;
	zz=0;
	for (int ii=0;ii<DIM*2;ii++)
		{
		if ((cc[ii]==0)&& (cc[ii+1]==0))
				buffer[zz]='A';
		else
			if ((cc[ii]==0)&&(cc[ii+1]==1))
				buffer[zz]='C';
			else
				if ((cc[ii]==1)&&(cc[ii+1]==0))
					buffer[zz]='G';
				else
					buffer[zz]='T';
			ii++;
			zz++;
	}
	buffer[zz]='\0';
};

virtual inline unsigned short getPoolCount(unsigned short pool) const {return 0;}
virtual inline unsigned short* getPools() const {return NULL;}
virtual inline unsigned short getPools(unsigned short pools[]) const {return 0;}
virtual inline unsigned short getNumPools() const {return 0;}
virtual inline unsigned short getNumBacs() const {return 0;}
virtual inline unsigned short getBacs(unsigned short bacs[]) const {return 0;}
//inline unsigned short getBac(unsigned short i) const {return 0;}
virtual void setBacs(const unsigned short bacs[], unsigned short numBacs){}

//!It compacts an unsigned integer in a char
void store_compact(unsigned int nval, char& c){
	register unsigned char cc;
	cc = (unsigned char)(0xFF & nval); 
	c=cc;
};

//!It uncompacts a char in an unsigned integer
void load_compact(char& c,unsigned int& nval)const{
	register unsigned char cc0 = c;
	register unsigned int uu = (unsigned int)(cc0 & 0xFF);
	nval=uu;
};


//!It encodes a bitset in an unsigned long int
unsigned long long int trans() const {
	register unsigned long long i=0;
	if (cc[(DIM*2)-1]==1)
		i++;
	for (int ii=(DIM*2)-2;ii>=0;ii--)
		{
		i=i<<1;
		if (cc[ii]==1)
			i++;
		}
	return i;
};

//!It encodes a bitset in a size_t
size_t trans2() const {
	register size_t i=0;
	if (cc[(DIM*2)-1]==1)
		i++;
	for (int ii=(DIM*2)-2;ii>=0;ii--)
		{
		i=i<<1;
		if (cc[ii]==1)
			i++;
		}
	return i;
}

void getBacs(unsigned short& bac1, unsigned short& bac2, unsigned short& bac3){};

};

class dataSet{
public:
	double v;
	int p;
	dataSet(){v=0.0; p=-1;};
	dataSet(double v,int p){this->v=v; this->p=p;};
//!Operator< 
	friend bool operator<(const dataSet& p1, const dataSet& p2){
  	
	if (p1.v<p2.v)
		return true;
	else		
		return false;
	};
};



class mybitsetBAC{
private:
	//!BAC ID.
	unsigned short id;
	//!Vectorbit for pools  
	char pools[POOLS2C];	
public:
	
//!Empty Constructor.
mybitsetBAC(){ 
	memset(pools,0, POOLS2C*sizeof(char) );
	id=MAXSHORT;
};
	
//!Copy Constructor mybitsetBAC
mybitsetBAC(mybitsetBAC& p){ 
	memset(pools,0, POOLS2C*sizeof(char) );
	for (int jj=0;jj<POOLS2C;jj++)
		{
		pools[jj]=p.pools[jj];
	}
	id=p.id;
};

//!Copy Constructor bitset
mybitsetBAC(bitset<POOLS>& p){
	 memset(pools,0, POOLS2C*sizeof(char) );
	for (int jj=0;jj<POOLS;jj++)
		{
		if (p[jj]==1)
			this->setPool(jj);
		}
	id=MAXSHORT;
};
	
//!It clears
void clear(){
	id=MAXSHORT;
	memset(pools, 0, POOLS2C*sizeof(char) );
};

inline void copy(mybitsetBAC& p) {};

//!Copy Constructor by mybitset
/*mybitsetBAC(mybitset& p){ 
	memset(pools,0, POOLS2C*sizeof(char) );
	for (int jj=0;jj<POOLS2C;jj++)
		{
		pools[jj]=p.getPool(jj);
		}
	id=MAXSHORT;
};*/

//!It sets to 1 the ith position  of the pools' bitvector.
inline void setPool(const int& pool){
	int ii=pool/8;
	int sh=pool%8;
	pools[ii]=pools[ii]|c[sh];
};

//!It resets to 0 the ith position  of the pools' bitvector.
	inline void resetPool(const int& pool){
	int ii=pool/8;
	int sh=pool%8;
	pools[ii]=pools[ii]&r[sh];
};

//!It returns the value of ith position  of  the pools' bitvector.
inline  unsigned int getPool(const int&i)const{
	int ii=i/8;
	int sh=i%8;
	if ((pools[ii] & c[sh]))
		return 1;
	else
		return 0;
};	

//!It sets to val the BAC id. 
inline void set(int val){
	id=val;
};

//!It returns the value of BAC id.
	inline int getId(){
	return id;
};


//!Return the number of 1 in the pool vectors
inline int count(void){
	int i=0;
	for (int jj=0;jj<POOLS2C;jj++)
		{
		for(int ii=0;ii<8;ii++)
			{
			if (pools[jj]&c[ii])
				{
				i++;
				}
			}
		}
	return i;
};


//!Operator==
friend bool operator==(mybitsetBAC& p1,mybitsetBAC& p2){
	for (int i=0;i<POOLS2C;i++)
		{
		if (p1.pools[i]!=p2.pools[i])
			{
			return false;
			}
		}
return true;
};


//!Operator<<
friend ostream& operator<<(ostream& out, const mybitsetBAC& p){
	p.print(out); 
	return out;
};

void print(ostream& out) const {
	//string buffer="";
	out<<id<<endl<<"%";
	for (int i=0;i<POOLS;i++) {
		if (getPool(i)==1)
			out<<" "<<i;
	}	
}

//!It encodes a bitset in an unsigned long int
unsigned long long int trans() {
	register unsigned int a=pools[0];
	for (int ii=1;ii<4;ii++)
		{
		a=a<<8;
		a+=pools[ii];
		}
	register unsigned int b=pools[4];
	for (int ii=5;ii<8;ii++)
		{
		b=b<<8;
		b+=pools[ii];
		}
	register unsigned int c=pools[8];
	for (int ii=9;ii<12;ii++)
		{
		c=c<<8;
		c+=pools[ii];
		}
 	a=a-b;  a=a-c;  a=a^(c >> 13);
  	b=b-c;  b=b-a;  b=b^(a << 8); 
  	c=c-a;  c=c-b;  c=c^(b >> 13);
  	a=a-b;  a=a-c;  a=a^(c >> 12);
  	b=b-c;  b=b-a;  b=b^(a << 16);
  	c=c-a;  c=c-b;  c=c^(b >> 5);
  	a=a-b;  a=a-c;  a=a^(c >> 3);
  	b=b-c;  b=b-a;  b=b^(a << 10);
  	c=c-a;  c=c-b;  c=c^(b >> 15);
	return (unsigned long long)c;
	};
	
};

};

