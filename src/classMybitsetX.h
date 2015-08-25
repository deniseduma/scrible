
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


/*#ifndef __BIT_H__
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
#endif*/

#ifndef __CMB_H__
	#define __CMB_H__
	#include "classMybitset.h"
#endif

#include <assert.h>

extern unsigned char c[8];
extern unsigned char r[8];
 
namespace MYBIT{

using namespace std;

class mybitsetx: public mybitset {
	
protected: 
//!Pool frequencies
unsigned short cw[POOLS]; 

public:


//!Empty Constructor.
mybitsetx(): mybitset() { 
	for (int i=0; i<POOLS; i++) 
		cw[i]=0;

};

//!Constructor
mybitsetx(const bitset<2*DIM>& b): mybitset(b) {
	//cw = new unsigned short[POOLS];
	for (int i=0; i<POOLS; i++) 
		cw[i]=0;
	//cw=(unsigned short*)malloc(POOLS*sizeof(unsigned short));
	//memset(cw,0, POOLS*sizeof(unsigned short));
}

//!Constructor
mybitsetx(const mybitset& b): mybitset(b) {
	for (int i=0; i<POOLS; i++) 
		cw[i]=0;
}

//!Copy constructor
mybitsetx(const mybitsetx& b): mybitset(b.cc) {
	//cw = new unsigned short[POOLS];
	for (int i=0; i<POOLS; i++) 
		cw[i]=b.cw[i];
	//cw=(unsigned short*)malloc(POOLS*sizeof(unsigned short));
	//memcpy(cw, b.cw, POOLS*sizeof(unsigned short));
}

//!Assignmenet operator
mybitsetx& operator=(const mybitsetx& b) {
	if (this!=&b) {
		cc=b.cc;
		for (int i=0; i<POOLS; i++) 
			cw[i]=b.cw[i];
	}
	return *this;
}

~mybitsetx() {
	//delete [] cw;
	//free(cw);
}
	
inline unsigned short getPool(unsigned int pool) const {
	return (cw[pool] > 0);
} 

inline unsigned short getPoolCount(unsigned short pool) const {
	return cw[pool];
}

inline unsigned short* getPools() {
	return cw;
}

//inline unsigned short getBac(unsigned short i) const {return 0;};
//inline unsigned short* getBacs() const {return NULL;};
//inline unsigned short getNumBacs() const {return 0;};
//inline void setBacs(const unsigned short bacs[], unsigned short numBacs){};

//!It increments cw by one.
inline void incPool(unsigned int pool){
	if(cw[pool]<MAXFREQ) 
		cw[pool]++;
};

//!Updates pool count by given value..
inline void updatePool(unsigned int pool, short val){
	if (((short)cw[pool] + val>=0)&&((short)cw[pool] + val<MAXFREQ))
		cw[pool]+=val;
};
	
//!It returns the value of cw.
inline unsigned short getPoolFreq(const int& pool) const{
	return cw[pool];
};	

//!It returns the number of pos pools.
inline unsigned short getNumPools() const {
	unsigned short numPools = 0;
	for (unsigned short pool=0; pool<POOLS; pool++)
		if (cw[pool]>0)
			numPools++;
	
	return numPools;		
};

//!It returns the pos pools and their number. 
inline unsigned short getPools(unsigned short pools[]) const {
	unsigned short numPools = 0;
	for (unsigned short pool=0; pool<POOLS; pool++)
		if (cw[pool]>0)
			pools[numPools++]=(pool+1);
	
	return numPools;		
};

//!Copy operator for pools' bitvector.
inline void copy(mybitsetx& p){
	for (unsigned int jj=0; jj<POOLS; jj++)
		p.cw[jj]=cw[jj];
};


//!Operator<<
/*friend ostream& operator<<(ostream& out, const mybitsetx& p){
	p.print(out);
	return out;
}*/

void print(ostream& out) const {
	
	mybitset::print(out); 
	out << "%" << endl;
	for (int i=0;i<POOLS;i++) {
		if (getPool(i)==1)
			out<<(i + 1)<<"("<<getPoolFreq(i)<<"), ";
	};
	out << endl;
};

};
};
