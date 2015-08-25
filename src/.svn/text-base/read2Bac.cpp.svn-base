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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <sys/time.h>
#include <sys/resource.h>


#ifndef __GNR_H__
        #define __GNR_H__
        #include "general.h"
#endif

#ifndef __CON_H__
        #define __CON_H__
       #include "const.h"
#endif

#ifndef __FSTREAM__
        #define __FSTREAM__
        #include <fstream>
#endif

#ifndef __SSTREAM__
        #define __SSTREAM__
        #include <sstream>
#endif

using namespace std;
using namespace general;

int main(int argc, char *argv[])
{
clock_t startGlobal,endGlobal;
startGlobal=clock();
ifstream in;

cout<<"\n\n =========================================================\n";
cout<<"|                   MAPPING READS IN BACS                  |\n";
cout<<" =========================================================\n";
cout<<"\n If you find any bug, send an email to beccuti@di.unito.it\n";

cout<<argc<<endl;
if (argc!=5)
        {
        std::cerr<<"\n\nUSE: read2BAC <path/ReadFile>  <path/OutputFile>  <number input files> <number output files>\n\n";
                exit(EXIT_FAILURE);
        }

int Ifiles=atoi(argv[3]);
int Ofiles=atoi(argv[4]);
string ReadFname(argv[1]);
string OFname=std::string(argv[2]);
clock_t startLocal;

cout<<"\n\n___________________________________________________________\n\n\t\t\tInput parameters \n___________________________________________________________\n\n";
cout<<"\tInput files name: "<<ReadFname<<"\n";
cout<<"\tOutput files name: "<<OFname<<"\n";
cout<<"\tInput files numbers: "<<Ifiles<<"\n";
cout<<"\tOutput files numbers: "<<Ofiles<<"\n";
cout<<"___________________________________________________________\n";



cout<<"\n\nSTART FORMATTING...\n"<<endl;

char buffer[MAXSIZE],name[MAXSIZE],bac[MAXSIZE],vec[MAXSIZE*2];
char delimC[] = " ";
Parser parser;
ofstream ** vecFp=(ofstream **) malloc ((Ofiles+1)*sizeof(ofstream*));
for (int i=0;i<Ofiles+1;i++)
        {
        ostringstream of;
        if (i>0)
                of<<OFname<<i;
        else
                of<<OFname<<"NoDec";
        vecFp[i]=new ofstream(of.str().c_str(),ofstream::out);
        if (!(*vecFp[i]))
                {
                cerr << "\n*****Error opening ouput file *****" << endl;
                exit(EXIT_FAILURE);
                }
        }
startGlobal=clock();
int num=0;
for (int i=0;i<Ifiles;i++)
	{
        int read=0;
        startLocal=clock();
        ostringstream of;
        of<<ReadFname<<i;
        cout<<"\nProcessing file: "<<of.str()<<endl;
        in.open(of.str().c_str(),ofstream::in);
        if(!in)
                {
                cerr << "\n*****Error opening input file *****" << endl;
                exit(EXIT_FAILURE);
                }
        while (!in.eof())
                {
                name[0]='\0';
                in.getline(name,MAXSIZE);
                num=in.gcount();
                if (name[num-1]!='\0')
                        {
                        name[num]='\0';
                        num++;
                        }
                if ((num>0)&&((name[0]=='@')||(name[0]=='>')))
                        {
                        buffer[0]='\0';
                        in.getline(buffer,MAXSIZE);
                        num=in.gcount();
                        if (buffer[num-1]!='\0')
                                {
                                buffer[num]='\0';
                                num++;
                                }
                        if (num==1)
                                {
                                cerr << "\n*****Error reading read "<<name<<" *****" << endl;
                                exit(EXIT_FAILURE);
                                }
                        in.getline(bac,MAXSIZE);
                        num=in.gcount();
                        if (bac[num-1]!='\0')
                                {
                                bac[num]='\0';
                                num++;
                                }
                        if (num==1)
                                {
                                cerr << "\n*****Error reading BAC of the read "<<name<<" *****" << endl;
                                exit(EXIT_FAILURE);
                                }
                        parser.update(delimC,bac);
                        read++;
                        if (read%1000000==0)
                                {
                                cout<<"Reads: "<<read<<endl;
                                }
                        bool dec=false;
                        for (unsigned int i=0; i< parser.size(); i++ )
                                {
                                int j=atoi(parser.get(i).c_str());
                                if (j!=-1)
                                        {
                                        if ((j<0||j>=Ofiles))
                                                {
                                                cerr << "\n*****Error BAC value for the read "<<name<<" *****" << endl;
                                                exit(EXIT_FAILURE);
                                                }
                                        else
                                                {
                                                *vecFp[j+1]<<name<<endl<<buffer<<endl;
                                                dec=true;
                                                }
                                        }
                                }
                        if (!dec)
                                {
                                *vecFp[0]<<name<<":$ ";
				in.getline(vec,2*MAXSIZE);
				num=in.gcount();
                                if (vec[num-1]!='\0')
                                        {
                                        vec[num]='\0';
                                        }

				if (vec[0]=='>')
					{
					cerr << "\n*****Error input file "<<name<<" *****" << endl;
                                        exit(EXIT_FAILURE);
					} 
				*vecFp[0]<<vec<<endl;
                                }
                        }

                }
        endGlobal=clock();
        cout<<"\n\tTime to read input file  "<<i<<":"<<((double)(endGlobal-startLocal))/CLOCKS_PER_SEC<<"s."<<endl;
        in.close();
        }

for (int i=0;i<Ofiles;i++)
        {
        (*vecFp[i]).close();
        delete(vecFp[i]);
	}
free(vecFp);

cout<<"\n\nEND FORMATTING\n"<<endl;
endGlobal=clock();
cout<<"\n=========================== TIME ===========================\n\n";
cout<<"\tTime to store reads in BAC files: "<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl;
cout<<"\n============================================================\n\n";

exit(EXIT_SUCCESS);
}
