#include<unistd.h> 
#include <stdio.h>
#include <vector>

#ifndef _CIO_H_
	#define _CIO_H_
	#include "classIO.h"
#endif

#define MAX_KMERS(x)  (100 - x + 1)

namespace IOFn {
using namespace general;
using namespace std;
using namespace MYBIT;

void readOverlaps(vector<char*>& overlaps, FILE* in);

/*class myhash//: public std::unary_function<mybitset, std::size_t> {
{
public:
	size_t operator()(const mybitsetx &b) const {
		return b.trans2();
	}
};


class myequal_to//: public std::binary_function<mybitset, mybitset, bool> {
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
*/

template <class A_Type>
void IOF<A_Type>::OverlapToKmers(int id, int l, int u) {

//FILE *in, *out;

clock_t startGlobal,endGlobal;

char delimC[] = ".";
Parser parser;
parser.update(delimC, ReadFname);

unsigned int num;

//char buffer[SIZEMAX];
char head[SIZEMAX];
char seq[SIZEMAX];


my_u_map kmers;
int* finalKmers = new int[POOLS*200]; 
int* kmersInPool = new int[POOLS];

for (int i=l;i<u;i++)
{
	ifstream in;
	ofstream out;
	ostringstream of, of1;
	
	of<<parser.get(0)<<i<<"."<<parser.get(1);
	//in=fopen(of.str().c_str(), "r");
	in.open(of.str().c_str(), ios_base::in);
	/*if (in==NULL) {
		cerr << "\nError opening input file " << of.str() << endl;
		exit(EXIT_FAILURE);
	}*/
	if(!in)
	{
		cerr<<"Error opening input file "<<of.str()<<endl;
		exit(EXIT_FAILURE);
	}

	of1<<OFname<<i<<".kmers" ;
	//out=fopen(of1.str().c_str(),"w");
	out.open(of1.str().c_str(), ios_base::out);
	/*if (out==NULL) {
		cerr << "\nError opening ouput file " << of1.str() << endl;
		exit(EXIT_FAILURE);
	}*/
	if(!out)
	{
		cerr<<"Error opening output file "<<of1.str()<<endl;
		exit(EXIT_FAILURE);
	}
	
	cout<<"["<<id<<"]"<< of.str()<< " "<<of1.str()<<endl;
	
	startGlobal=clock();
	unsigned int lowKmersReads = 0;
	
	//while (!feof(in)) {
	while (!in.eof()) {
		
		kmers.clear();
		//finalKmers.clear();
		for (int j=0; j<POOLS; j++)
			kmersInPool[j] = 0;

		head[0]='\0';
		//if (fgets(head, SIZEMAX, in)==NULL) 
		//	break;
		in.getline(head, SIZEMAX);
		
		num = in.gcount();
		//num=strlen(head);
		/*if (head[num-1]=='\n')
		{
			head[num-1]='\0';
			num--;
		}*/
		num--;
		
		if ((num>0)&&((head[0]=='@')||(head[0]=='>'))) {
		
		vector<string> reads;
		
		string token;
		istringstream is (head);
		while (getline(is, token, ' ')) 
			reads.push_back(token);
		
		for (unsigned int r=0; r<reads.size(); r++) {
			seq[0]='\0';
			//if (fgets(seq, SIZEMAX, in)==NULL) 
			//	break;
			in.getline(seq, SIZEMAX);
			num = in.gcount();
			//num=strlen(seq);
			/*if (seq[num-1]=='\n')
			{
				seq[num-1]='\0';
				num--;
			}*/
			num--;
			
			assert((num>=75)&&((seq[0]=='A')||(seq[0]=='G')||(seq[0]=='C')||(seq[0]=='T'))); 
			//if ((num>=kmerSize)&&((seq[0]=='A')||(seq[0]=='G')||(seq[0]=='C')||(seq[0]=='T')||(seq[0]=='N'))) {

			class mybitsetx b,d,e;
			for (unsigned int k=0; k<=num-kmerSize+1; k++) {
				int z=0;
				int j=kmerSize*2-1;
				bool findN=false;
				for (unsigned int c=0; c<kmerSize; c++) {
              				switch (seq[c+k]) {
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
					unsigned int count  = 1;
					auto it = kmers.find(e);
					if (it != kmers.end()) 
						count += it->second;
					kmers[e]=count;
				}	
			}//end for k
		//}//end if	
	}//end for r	
	
	/*sort kmers in reverse order of their counts within the cluster*/
	my_m_map counts;
	for (auto& it: kmers)
		counts.insert(pair<unsigned int, mybitsetx>(it.second, it.first));
	
	//DEBUG
	//if (kmers.size()>=500)
	//	printf("\n\t%s, num reads %lu, distinct kmers %lu\n", reads[0], reads.size(), kmers.size());
	unsigned int numKmers =0;
	for (auto& it: counts) {
		mybitsetx b = it.second;
		unsigned int count = it.first;
		
	if ((reads.size()==1) || (reads.size()==2 && count==2) || (reads.size()==3 && count>=2) || (count >= 3)) 
	{ 
		if (sh.searchCopy(b)) {
				
			for (int pool = 0; pool<POOLS; pool++) {
				finalKmers[numKmers*POOLS + pool] = b.getPoolFreq(pool);
				kmersInPool[pool] += b.getPool(pool);
			}
			
			numKmers++;
			if (numKmers==75)
				break;
		}
	}	
	}//end for
		
	for (unsigned int r=0; r<reads.size(); r++) 
		//fprintf(out, "%s ", reads[r]);
		out << reads[r] << " ";
	//fprintf(out, "kmers: %d bacs: \n", numKmers);
	out << "kmers: " << numKmers << " bacs: " <<  endl;

	if (numKmers<3) {
		
		//DEBUG
		/*printf("reads size %lu, total kmers %lu, numKmers %d\n", 
			reads.size(), kmers.size(), numKmers);
		for (unsigned int r=0; r<reads.size(); r++) 
			printf("%s ", reads[r]);
		printf("\n");
		if (reads.size()==1) 
			printf("%s\n", seq);
		if (kmers.size() < 3) { 
			for (auto& it2 : kmers) 
				cout << it2.first << it2.second << endl;
			printf("\n");
		}*/
		
		lowKmersReads++;
		
		continue;
	}

	for (unsigned int pool=0; pool<POOLS; pool++) {
		if (kmersInPool[pool]>0) {
			//fprintf(out, "%d ", pool);
			out<<pool<<" ";
			for (unsigned int k=0; k<numKmers; k++) {
				//fprintf(out, "%d ", finalKmers[k*POOLS + pool]);
				out<<finalKmers[k*POOLS + pool]<<" ";
			}	
			//fprintf(out, "\n");
			out<<endl;
		}
	}
	
	}//if ((num>0)&&((head[0]=='@')||(head[0]=='>'))) {
	
	}//end while (!feof(in))


endGlobal=clock();
cout << "\n[overlap2kmers:Overlap2Kmers]Time to process pool "<<i<<": "<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl;
cout << "Number of reads with less than 3 distinct kmers is " << lowKmersReads << endl << endl;

//fclose(in);
//fclose(out);
in.close();
out.close();

}//end for (int i=l;i<=u;i++)

delete [] finalKmers;
delete [] kmersInPool;

};


template <class A_Type>
void IOF<A_Type>::OverlapToKmersLight(int id, int l, int u) {

//FILE *in, *out;

clock_t startGlobal,endGlobal;

char delimC[] = ".";
Parser parser;
parser.update(delimC, ReadFname);

unsigned int num;

char head[SIZEMAX];
char seq[SIZEMAX];

my_u_map kmers;
for (int i=l;i<=u;i++)
{
	ifstream in;
	ofstream out;
	ostringstream of, of1;
		
	of<<parser.get(0)<<i<<"."<<parser.get(1);
	//in=fopen(of.str().c_str(), "r");
	in.open(of.str().c_str(), ios_base::in);
	/*if (in==NULL) {
		cerr << "\nError opening input file " << of.str() << endl;
		exit(EXIT_FAILURE);
	}*/
	if(!in)
	{
		cerr << "Error opening input file " << of.str().c_str() << endl;
		exit(EXIT_FAILURE);
	}

	of1<<OFname<<i;
	//out=fopen(of1.str().c_str(),"w");
	out.open(of1.str().c_str(),ios_base::out);
	/*if (out==NULL) {
		cerr << "\nError opening ouput file " << of1.str() << endl;
		exit(EXIT_FAILURE);
	}*/
	if (!out)
	{
		cerr << "Error opening output file " << of1.str().c_str() << endl;
		exit(EXIT_FAILURE);
	}
	
	cout << "[" << id << "]" << of.str().c_str() << " " << of1.str().c_str() << endl;

	startGlobal=clock();
	
	unsigned int splitClusters = 0;
	
	//while (!feof(in)) {
	while (!in.eof()) {
		
		kmers.clear();

		head[0]='\0';
		//if (fgets(head, SIZEMAX, in)==NULL) 
		//	break;
		in.getline(head,SIZEMAX);
		
		//num=strlen(head);
		num=in.gcount();
		/*if (head[num-1]=='\n')
		{
			head[num-1]='\0';
			num--;
		}*/ 
		num--;

		if ((num>0)&&((head[0]=='@')||(head[0]=='>'))) {
		
		//DEBUG
		//cout << "head:" << head << endl;
		vector<string> reads;
		vector<string> seqs;
		
		string token;
		istringstream is (head);
		while (getline(is, token, ' ')) { 
			//cout << token << endl;
			reads.push_back(token);
		}	
	
		//DEBUG
		/*cout << "reads: ";
		for (unsigned int r=0; r<reads.size(); r++) 
			cout << reads[r] << " ";
		cout << endl << endl;*/	

		for (unsigned int r=0; r<reads.size(); r++) {
			seq[0]='\0';
			//if (fgets(seq, SIZEMAX, in)==NULL) 
			//	break;
			in.getline(seq,SIZEMAX);
			
			//num=strlen(seq);
			num=in.gcount();
			/*if (seq[num-1]=='\n')
			{
				seq[num-1]='\0';
				num--;
			}*/ 
			num--;
			seqs.push_back(seq);
			
			assert ((num>=75)&&((seq[0]=='A')||(seq[0]=='G')||(seq[0]=='C')||(seq[0]=='T'))); 
			//if ((num>=kmerSize)&&((seq[0]=='A')||(seq[0]=='G')||(seq[0]=='C')||(seq[0]=='T')/*||(seq[0]=='N')*/)) {

			class mybitsetx b,d,e;
			for (unsigned int k=0; k<=num-kmerSize+1; k++) {
				int z=0;
				int j=kmerSize*2-1;
				bool findN=false;
				for (unsigned int c=0; c<kmerSize; c++) {
              				switch (seq[c+k]) {
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
					unsigned short count  = 1;
					auto it = kmers.find(e);
					if (it != kmers.end()) 
						count += it->second;
					kmers[e]=count;
				}	
			}//end for k
		//}//end if 
		/*else {
			cout<<"assertion broken"<<endl; 
			cout<<head<<endl;
			cout<<"the sequence"<<endl;
			cout<<seq<<endl<<endl;
			exit(1);
		}*/
	}//end for r	
	
	unsigned int inHash =0;
	unsigned int numKmers =0;
	for (auto& it: kmers) {
		mybitsetx b = it.first;
		unsigned int count = it.second;
		
		if ((reads.size()==1) || (reads.size()==2 && count==2) || (reads.size()==3 && count>=2) || (count >= 3)) { 
			numKmers++; 
			if (sh.search(b))
				inHash++;
		}
	}
        
	if ((reads.size()==1) || (inHash>=3)) {
		for (unsigned int r=0; r<(reads.size() - 1); r++) 
			//fprintf(out, "%s ", reads[r]);
			out << reads[r] << " ";
		//fprintf(out, "%s\n", reads[reads.size() - 1]);	
		out << reads[reads.size()-1] << endl;
	} else {//split bad clusters
		//fprintf(out, "%s\n", reads[0]);
		out << reads[0] << endl;
		for (unsigned int r=1; r<reads.size(); r++) 
			//fprintf(out, ">%s\n", reads[r]);
			out << ">" << reads[r] << endl;
		
		//DEBUG
		printf("Split cluster of size %lu, having kmers size %lu, numKmers %d and inHash %d\n", reads.size(), kmers.size(), numKmers, inHash);	
		cout << head << endl;
		cout << "reads" << endl;
		for (unsigned int r=0; r<reads.size(); r++) 
			cout << reads[r] << endl << seqs[r] << endl;
		//cout << "kmers" << endl;
		//for (auto& it : kmers)	
		//	cout<<it.first<<cout<<it.second<<endl;
		
		splitClusters++;	
	}
	}//end if if ((num>0)&&((head[0]=='@')||(head[0]=='>')))	
	}//end for o


endGlobal=clock();
printf("\n[overlap2kmers:OverlapToKmersLight]Time to process pool %d: %.2fs \n[%d]Number of bad clusters which were split is %d\n\n", i, ((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC, id, splitClusters);

//fclose(in);
//fclose(out);
in.close();
out.close();

}//end for (int i=l;i<=u;i++)

};


template <class A_Type>
void IOF<A_Type>::RemoveLowComplexity(int id, int l, int u) {

FILE *in, *out, *err;

clock_t startGlobal,endGlobal;

char delimC[] = ".";
Parser parser;
parser.update(delimC, ReadFname);

unsigned int num;

char head[SIZEMAX];
char seq[SIZEMAX];

my_u_map kmers;

for (int i=l;i<=u;i++)
{
	ostringstream of, of1, of2;
	
	of<<parser.get(0)<<i<<"."<<parser.get(1);
	in=fopen(of.str().c_str(), "r");
	if (in==NULL) {
		cerr << "\nError opening input file " << of.str() << endl;
		exit(EXIT_FAILURE);
	}

	of1<<OFname<<i<<".complex" ;
	out=fopen(of1.str().c_str(),"w");
	if (out==NULL) {
		cerr << "\nError opening ouput file " << of1.str() << endl;
		exit(EXIT_FAILURE);
	}
	
	of2<<OFname<<i<<".lowcomplex" ;
	err=fopen(of2.str().c_str(),"w");
	if (err==NULL) {
		cerr << "\nError opening discarded reads file " << of1.str() << endl;
		exit(EXIT_FAILURE);
	}
	
	unsigned int lowComplexityReads = 0;
	
	startGlobal=clock();
	
	while (!feof(in)) {
		
		kmers.clear();

		head[0]='\0';
		if (fgets(head, SIZEMAX, in)==NULL) 
			break;
		
		num=strlen(head);
		if (head[num-1]=='\n')
		{
			head[num-1]='\0';
			num--;
		}
		
		//char buf[num + 1];
		//strncpy(buf, head, num);

		vector<char*> reads;
		
		char* item = strtok(head, " \n");
		while (item != NULL) {
			reads.push_back(item);
			item = strtok(NULL, " \n");
		}

		for (unsigned int r=0; r<reads.size(); r++) {
			seq[0]='\0';
			if (fgets(seq, SIZEMAX, in)==NULL)
				break;
			num=strlen(seq);
			if (seq[num-1]=='\n')
			{
				seq[num-1]='\0';
				num--;
			}
						
			if ((num>=kmerSize)&&((seq[0]=='A')||(seq[0]=='G')||(seq[0]=='C')||(seq[0]=='T')||(seq[0]=='N'))) 
			{

			class mybitsetx b,d,e;
			for (unsigned int k=0; k<=num-kmerSize+1; k++) {
				int z=0;
				int j=kmerSize*2-1;
				bool findN=false;
				for (unsigned int c=0; c<kmerSize; c++) {
              				switch (seq[c+k]) {
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
					unsigned short count  = 1;
					auto it = kmers.find(e);
					if (it != kmers.end()) 
						count += it->second;
					kmers[e]=count;
				}	
			}//end for k
		}//end if	
	}//end for r	
	
	//DEBUG
	/*if ((!strcmp(reads[0], ">a0001A06_5/1")) || 
		(!strcmp(reads[0], ">a0001A06_12ef/1")) ||
		(!strcmp(reads[0], ">a0001A06_31/1")) || 
		(!strcmp(reads[0], ">a0001A06_e10/1")) ||
		(!strcmp(reads[0], ">a0001A06_e11/1")) ||
		(!strcmp(reads[0], ">a0001A06_fdf/2"))) {
		printf("%s\n%s\ndistinct %lu,total %d\n", reads[0], seq, 
			kmers.size(), num-kmerSize+1);
	}*/
	
	if (kmers.size()/((double)(num-kmerSize+1)) <= 0.5) {
		for (unsigned int r=0; r<reads.size(); r++) 
			fprintf(err, "%s\n%s\n", reads[r], seq);
		lowComplexityReads++;
	} else {
		for (unsigned int r=0; r<reads.size(); r++) 
			fprintf(out, "%s\n%s\n", reads[r], seq);
	}
	
	}//end for o


endGlobal=clock();
cout << "\n[overlap2kmers:RemoveLowComplexity]Time to process pool "<<i<<": "<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl;
cout << "Number of discarded low complexity reads is " << lowComplexityReads << endl << endl;

fclose(in);
fclose(out);
fclose(err);

}//end for (int i=l;i<=u;i++)

};


void readOverlaps(vector<char*>& overlaps, FILE* in) 
{
int num;
int count = 0;
char* buffer;	
//char buffer[SIZEMAX];
while (!feof(in)) {
	buffer=new char[SIZEMAX];	
	if (fgets(buffer, SIZEMAX, in)==NULL) {
		break;
	}	
	num=strlen(buffer);
	if (buffer[num-1]=='\n')
	{
		buffer[num-1]='\0';
		num--;
	}
	overlaps.push_back(buffer);
	
	if (buffer[0]=='>') 
		count++;

}
	//DEBUG
	printf("[readOverlaps]num overlaps read from file %d\n", count);
};


template <class A_Type>
int IOF<A_Type>::readOverlaps(int i, struct overlap* overlaps, ofstream& out) {

ostringstream of, of1;

char delimC[] = ".";
Parser parser;
parser.update(delimC, ReadFname);

ifstream in;
of<<parser.get(0)<<i<<"."<<parser.get(1);
//in=fopen(of.str().c_str(), "r");
in.open(of.str().c_str(), ios_base::in);
if(!in) {
	cerr<<"Error opening input file "<<of.str()<<endl;
	exit(EXIT_FAILURE);
}

parser.update(delimC, OFname);

of1<<parser.get(0)<<i<<"."<<parser.get(1)<<"."<<parser.get(2);
out.open(of1.str().c_str(), ios_base::out);
if(!out)
{
	cerr<<"Error opening output file "<<of1.str()<<endl;
	exit(1);
}

//DEBUG
printf("Reading overlaps from %s\n", of.str().c_str());
printf("Output file is %s\n", of1.str().c_str());

unsigned int num;
char head[SIZEMAX];
char seq[READSIZE+1];

int index=-1;
int numOverlaps=0;
struct overlap* crt;

clock_t startGlobal, endGlobal;

startGlobal=clock();
//while (!feof(in)) {
while (in.good()) {
	
	head[0]='\0';
	in.getline(head, SIZEMAX);
	num=strlen(head);
	
	//DEBUG
	//if (num>=1000000) {
	//	printf("%s\n", head);
	//	printf("length of header %d\n\n", num);
	//}	
	
	if ((num>0)&&((head[0]=='@')||(head[0]=='>'))) {
		
		index++;
		crt = overlaps + index;
		
		numOverlaps++;
	
		crt->header=(char*)malloc((num+1+1000)*sizeof(char));
		strncpy(crt->header, head, num);
		crt->header[num]='\0';
	
		//DEBUG
		//printf("%s\n", crt->header);
		//printf("numOverlaps so far %d\n", numOverlaps);
	
		string token;
		vector<string> reads;
		istringstream is (head);
		while (getline(is, token, ' ')) 
			reads.push_back(token);
		
		crt->numReads=reads.size();
		crt->reads=(char**)malloc(((crt->numReads)*sizeof(char*)));
	
		for (int r=0;r<(crt->numReads);r++) {
			seq[0]='\0';
			in.getline(seq, READSIZE+1);
			num=strlen(seq);
					
			assert((num>=75)&&((seq[0]=='A')||(seq[0]=='G')||(seq[0]=='C')||(seq[0]=='T'))); 
		
			crt->reads[r]=(char*)malloc((num+1)*sizeof(char));
			strncpy(crt->reads[r], seq, num);
			crt->reads[r][num]='\0';
			
			//DEBUG
			//printf("%s\n", crt->reads[r]);
	
		}//end for r
	}//end if
}//end while

in.close();

endGlobal=clock();
cout << "\n[overlap2kmers:ReadOverlaps]Time to process pool "<<i<<": "<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl;
	
return numOverlaps;

};


template <class A_Type>
void IOF<A_Type>::OverlapToDeconv(int id, int readNo, const SpMatCol& phi, int* bacPools, 
	struct overlap* overlap, ofstream& out)
{

using namespace Eigen;

my_u_map kmers;

SpMatCol finalKmers(POOLS, MAX_KMERS(kmerSize));
vector<TripletF> tripletList(POOLS*MAX_KMERS(kmerSize));
//tripletList.reserve(POOLS*MAX_KMERS(kmerSize));

char* header=overlap->header;
int numReads=overlap->numReads;
char** reads=overlap->reads;

//DEBUG
//if (readNo%100000==0)
//	printf("Thread %d read number %d\n", id, readNo);
//	printf("%s\n\n", header);

//clock_t startGlobal, endGlobal;

//startGlobal=clock();
for (int r=0;r<numReads;r++) {
	char* read = reads[r];
	int num = strlen(read);
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
			unsigned int count  = 1;
			auto it = kmers.find(e);
			if (it != kmers.end()) 
				count += it->second;
				kmers[e]=count;
		}	
	}//end for k
}//end for r	
	
my_m_map counts;
for (auto& it: kmers)
	counts.insert(pair<unsigned int, mybitsetx>(it.second, it.first));
	

unsigned int numKmers =0;
for (auto& it: counts) {
	mybitsetx b = it.second;
	unsigned int count = it.first;
		
	if ((numReads==1)||(numReads==2&&count==2)||(numReads==3&&count>=2)||(count>=3)){ 
		if (sh.searchCopy(b)) {
				
			for (int pool = 0; pool<POOLS; pool++) //{
				//finalKmers[numKmers*POOLS+pool]=b.getPoolFreq(pool);
				//kmersInPool[pool]+=b.getPool(pool);
				tripletList.push_back(TripletF(pool,numKmers,b.getPoolFreq(pool)));
			//}
			
			numKmers++;
			if (numKmers==75)
				break;
		}
	}	
}//end for

finalKmers.setFromTriplets(tripletList.begin(), tripletList.end());

//deconv the overlap
//deconvOverlap(phi, bacPools, header, numKmers, finalKmers, out);

//DEBUG
//endGlobal=clock();
//if (readNo%100000==0)
//	cout << "\n[overlap2kmers:OverlapToDeconv]Time to process read "<<readNo<<": "<<((double)(endGlobal-startGlobal))/(CLOCKS_PER_SEC)<<"s."<<endl;

};
}

