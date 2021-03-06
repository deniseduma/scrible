/***************************************************************************
 *   Copyright (C) 2011 by Marco Beccuti   *
virtual ~mysig() {}
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

#include <string.h>
#include <stdlib.h>

extern unsigned char c[8];
extern unsigned char r[8];
extern unsigned int hashSizes[SIZES];

namespace MYBIT{

using namespace std;

class mysig {

protected:
	//!A sig (set of pools) and corresponding bacs.
	unsigned short* bacs;
	unsigned short* pools;
	unsigned short numPools, numBacs; 
public:

//!Empty Constructor.
mysig(): bacs(new unsigned short[MAX_BACS]), pools(new unsigned short[POOLS]), numPools(0), numBacs(0) {};

//!Copy constructor
mysig(const mysig& s) {
	numBacs=s.numBacs;
	numPools=s.numPools;
	bacs=new unsigned short[numBacs];
	for (unsigned short i=0;i<numBacs;i++)
		bacs[i]=s.bacs[i];
	pools=new unsigned short[numPools]; 
	for (unsigned short i=0;i<numPools;i++)
		pools[i]=s.pools[i];
}

~mysig() {
	delete [] bacs; 
	bacs=NULL; 
	delete [] pools; 
	pools=NULL;
}

//!Assignmenet operator
mysig& operator=(const mysig& s) {
	if (this!=&s) {
		unsigned short i;
		numPools=0; numBacs=0;
		delete bacs; delete pools;
		bacs=new unsigned short[s.numBacs];
		pools=new unsigned short[s.numPools];
		for (i=0;i<s.numPools;i++)
			pools[numPools++]=s.pools[i];
		for (i=0;i<s.numBacs;i++)
			bacs[numBacs++]=s.bacs[i];
	}	
	return *this;
}
	
//!It returns the bacs
inline unsigned short* getBacs() const {
	return bacs;
};

inline unsigned short* getPools() const {
	return pools;
};

inline unsigned short getNumPools() const {
	return numPools;
}

inline unsigned short& getNumPools() {
	return numPools;
}

inline unsigned short getNumBacs() const {
	return numBacs;
}

inline unsigned short& getNumBacs() {
	return numBacs;
}

inline unsigned short getBac(unsigned int i) const {
	return bacs[i];
}

inline unsigned short getPool(unsigned int i) const {
	return pools[i];
}

inline void resetNumBacs() {
	numBacs=0;
}

inline void setPools(unsigned short* pools, unsigned short numPools) {
	this->numPools=numPools;
	for (unsigned short i=0;i<numPools;i++)
		this->pools[i]=pools[i];
}

inline void setBacs(unsigned short* bacs, unsigned short numBacs) {
	this->numBacs=numBacs;
	for (unsigned short i=0;i<numBacs;i++)
		this->bacs[i]=bacs[i];
}

inline void copy(mysig& s) const {
	s.setBacs(bacs, numBacs);
};


//!It compares the two sequences  p1 and p2  and returns 0 when p1=p2, 1 when p1>p2 and 2 when p1<p2. 
/*inline int compare(const mysig& s2){
	unsigned int i=0;
 	while (((this->pools[i]==s2.pools[i]))&&(i<numPools))
 			{
 			i++;		
 			}
 	if (i==numPools)
 			{
 			return 0;
 			}
 	if (this->cc[i]==1)
 			return 1;
 		else		
 			return 2;
};*/


//!Operator==
friend bool operator==(const mysig& s1, const mysig& s2){
	if (s1.numPools!=s2.numPools)
		return false;
 	
	unsigned short i=0;
	while (((s1.pools[i]==s2.pools[i]))&&(i<s1.numPools))
 		i++;		
 	
	if (i==s1.numPools)
		return true;
	
	return false;
};

//!Operator!=
/*friend bool operator!=(mybitset& p1,mybitset& p2){
	if (p1.cc!=p2.cc) {
		return true;
	};
	return false;
};*/


//!Operator< 
/*friend bool operator<(const mybitset& p1, const mybitset& p2){
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
};*/


//!Operator<<
friend ostream& operator<<(ostream& out, const mysig& s){
	s.print(out);
	return out;
}

virtual void print(ostream& out) const { 
	out << "sig" << endl;
	out << "numPools " << numPools << endl;
	out << "pools: ";
	unsigned short i;
	for (i=0;i<numPools;i++)
		out << pools[i] << " ";
	out << endl;
	out << "numBacs " << numBacs << endl;
	out << "bacs: ";
	for (i=0;i<numBacs;i++)
		out << bacs[i] << " ";
};


//!It encodes a bitset in an unsigned long int
unsigned long long trans() const {
	register unsigned long long t=pools[0];
	for (unsigned short i=1;i<numPools;i++) {
		//if (pools[i]%2==1)
		//	t=t<<1;
		t*=100;
		t+=pools[i];
	}	
	return ((A*t+L)%hashSizes[2]);
};

//!It encodes a bitset in a size_t
/*size_t trans2() const {
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
}*/

};

};

