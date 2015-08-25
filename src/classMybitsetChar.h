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

extern unsigned char c[8];
extern unsigned char r[8];


 
namespace MYBIT{

using namespace std;



class mybitsetChar{
//!Sequence  in binary format.
bitset <DIM*2> cc;
//!Vectorbit for pools  
char pp[POOLS2C];
#if REPET
//unsigned short cw[NUMC]; //barley
unsigned char  cw[NUMC]; //rice

#endif	
public:
	unsigned int Inext;
	//!Empty Constructor.
	mybitsetChar(){ 
	memset(pp,0, POOLS2C*sizeof(char) );
#if REPET
	//memset(cw,0, NUMC*sizeof(unsigned short) );//barley
	memset(cw,0, NUMC*sizeof(unsigned char) ); //rice
	cw[0]=1;
#endif		
	};
	
	inline bitset<2*DIM> getBitset() {
		return cc;
	}

	//!It sets to 1 the ith position  of the pools' bitvector.
	inline void setVP(const int& i){
	int ii=i/8;
	int sh=i%8;
	pp[ii]=pp[ii]|c[sh];
	};
	//!It resets to 0 the ith position  of the pools' bitvector.
	inline void resetVP(const int& i){
	int ii=i/8;
	int sh=i%8;
	pp[ii]=pp[ii]&r[sh];
	};
	//!It returns the value of ith position  of  the pools' bitvector.
	inline unsigned int getVP(const int&i)const{
	int ii=i/8;
	int sh=i%8;
	if ((pp[ii] & c[sh]))
		return 1;
	else
		return 0;
	};
	
	
	//!It returns the value of ith position of pools vector in char format (low-level)
	inline char getVPChar(const int&i) {
	if ((i<0)|| (i>POOLS2C))
		return '\0';
	return pp[i];
	};

	//!It sets to val  the ith position  of the sequence's bitvector. 
	inline void set(const int& i,int val){
	cc[i]=val;
	};
	//!It returns the value of ith position of the sequence's bitvector.
	inline int get(const int& i){
	return cc[i];
	};

#if REPET	
	//!It increments cw by one.
	inline void inc(){

	int p1bit=this->count();
	if((p1bit>0)&&(p1bit<=NUMC)&&(this->cw[p1bit-1]<MAXFREQ))
		{
		this->cw[p1bit-1]++;
		}
	};
	
#endif

#if REPET
	//!It returns the value of cw.
	inline unsigned int getCW(int i)const{
	return cw[i];
	};	
#endif	
	//!Copy operator for pools' bitvector.
	inline int copy(mybitsetChar& p){
	int i=0;
	for (int jj=0;jj<POOLS2C;jj++)
		{
		p.pp[jj]=pp[jj];
		for(int ii=0;ii<8;ii++)
			{
			if (pp[jj]&c[ii])
				{
				i++;
				}
			}
		}
#if REPET
	for (int jj=0;jj<NUMC;jj++)
		p.cw[jj]=cw[jj];
#endif
	return i;
	};

	//!Return the number of 1 in the pool vectors
	inline int count(void){
	int i=0;
	for (int jj=0;jj<POOLS2C;jj++)
		{
		for(int ii=0;ii<8;ii++)
			{
			if (pp[jj]&c[ii])
				{
				i++;
				}
			}
		}
	return i;
	};

	//!Return the number of 1 in the union of pools vectors
        inline int count(const mybitsetChar& p){
        int i=0;
        for (int jj=0;jj<POOLS2C;jj++)
                {
                for(int ii=0;ii<8;ii++)
                        {
                        if ((pp[jj]&c[ii])||(p.pp[jj]&c[ii]))
                                {
                                i++;
                                }
                        }
                }
        return i;
        };

	//!Check if they share at least one pool in their pool bitvector. 
	inline  bool EmptyIntersection(const mybitsetChar& p){
	for (int jj=0;jj<POOLS2C;jj++)
		{
		for(int ii=0;ii<8;ii++)
			{
			if ((pp[jj]&c[ii])&&(p.pp[jj]&c[ii]))
				{
				return false;
				}
			}
		}
	return true;
	};
	//!Union operator between the two pools' bitvector. The results is stored in p.
	inline void unionPV(bitset <POOLS>& p){
	register char v;
	unsigned int kk=0;
	for (int jj=0;jj<POOLS2C;jj++)
		{
		v=pp[jj];
		for(int ii=0;ii<8;ii++)
			{
			if (v&c[ii])
				{
				p[kk]=1;
				}
			kk++;
			}
		}
	};
	//!Union operator between the two pools' bitvector encoding the pools frequency on vector fv. The results is stored in p.
	inline void unionPV(bitset <POOLS>& p,unsigned int fv[]){
	register char v;
	int kk=0;
	for (int jj=0;jj<POOLS2C;jj++)
		{
		v=pp[jj];
		for(int ii=0;ii<8;ii++)
			{
			if (v&c[ii])
				{
				p[kk]=1;
				fv[kk]++;
				}
			kk++;
			}
		}
	};

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
	}
	//!It compares the two sequences  p1 and p2  and returns 0 when p1=p2, 1 when p1>p2 and 2 when p1<p2. 
	inline int compare(const mybitsetChar& p2){
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
	}
	//!Operator==
	friend bool operator==(mybitsetChar& p1,mybitsetChar& p2){
	if (p1.cc==p2.cc)
		{
		return true;
		}
	else
		return false;
	};


	//!Operator< 
	friend bool operator<(const mybitsetChar& p1, const mybitsetChar& p2){
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
	friend ostream& operator<<(ostream& out, const mybitsetChar& p){
	int zz=0;
	string buffer="";
	zz=0;
	for (int ii=0;ii<DIM*2;ii++)
		{
		if ((p.cc[ii]==0)&& (p.cc[ii+1]==0))
				buffer+='A';
		else
			if ((p.cc[ii]==0)&&(p.cc[ii+1]==1))
				buffer+='C';
			else
				if ((p.cc[ii]==1)&&(p.cc[ii+1]==0))
					buffer+='G';
				else
					buffer+='T';
			ii++;
			zz++;
		}
	
	out<<buffer<<endl<<"%";
	for (int i=0;i<POOLS;i++)
		{
		if (p.getVP(i)==1)
			out<<" "<<i;
		//out<<p.getVP(i);
		}
#if DEBUG 
	if (p.Inext==DEFAULTP)
		out<<"\nNext: NULL"<<endl;
	else
		out<<"\nNext: "<<p.Inext<<endl;
	
#endif
	return out;
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
	unsigned long long int trans()
	{
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
	
};

	

}

