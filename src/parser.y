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

%{
#ifndef __OPER_H__
	#define __OPER_H__
	#include "operation.h"
#endif



#define DEBUGPARSER 0

extern "C"
{
#include "parser.h"
#include <stdio.h>
#include <string.h>

	int yylex(void);
        int yywrap()
        	{
                return 1;
        	}
	void yyerror(const char *str)
		{
        	fprintf(stderr,"errore: %s\n",str);
		}
}
extern FILE * yyin; //standard yacc from file
int yyparse();

using namespace OPER;


class OPERATOR tmp;
int type;
int i=0;
int setop=0;
using namespace std;
vector <class OPERATOR>* oper;

int initParser(std::string filename,vector <class OPERATOR>& operations)
{
oper=&operations;
std::string file=filename+std::string(".operation");
yyin = fopen(file.c_str(),"r");
if (yyin==NULL)
    {
    printf("\n\nError: opening input file: %s\n\n",file.c_str());
    exit(EXIT_FAILURE);
    }
yyparse();
fclose(yyin);
#if DEBUGPARSER
vector<class OPERATOR>::iterator iterm = operations.begin();
while (iterm!=operations.end())
	{
	cout<<"****************************************\n";
	cout<<"\t\tType: "<<iterm->type<<endl;
	vector <int>::iterator iterl = iterm->p1.begin();
	while (iterl!=iterm->p1.end())
		{
		cout<<"\t Op1s: "<<(*iterl)<<" ";
		iterl++;
		}
	cout<<endl;
	iterl = iterm->p2.begin();
	while (iterl!=iterm->p2.end())
		{
		cout<<"\t Op2s: "<<(*iterl)<<" ";
		iterl++;
		}
	cout<<endl;	
	iterm++;
	}
	cout<<endl;
#endif
return EXIT_SUCCESS;
}

%}






//YACC GRAMMAR
%union{
int num;
char var[255];
}

%token <num> NUMBER
%token <var> STRING TO AVERAGE CORREC KMERCOUNT

%%
File:	File Formula ';'
	|
	Formula ';'
	;

Formula: 
	Av  
	{
	tmp.type = type;
	oper->push_back(tmp);
	tmp.p1.clear();
	tmp.p2.clear();
	}
	|
	Co 
	{
	tmp.type = type;
	oper->push_back(tmp);
	tmp.p1.clear();
	tmp.p2.clear();
	}
	|
	KmC	
	{
	tmp.type = type;
	oper->push_back(tmp);
	tmp.p1.clear();
	tmp.p2.clear();
	}
	;

Av:    AVERAGE List 
		{
		type = AVERAGE;
		setop=0;
		}
	;

Co:
	CORREC List To  
		{
		type= CORREC;
		setop=0;
		} 
	;

KmC:	KMERCOUNT
		{
		type= KMERCOUNT;
		setop=0;
		} 
	;
To:	STRING List
	;

List:   List ',' Op
		{
		setop=1;
		}
	|
	Op
		{
		setop=1;
		}
	;


Op: NUMBER 
		{
		
		if (setop)
			tmp.p2.push_back($1);
		else
			tmp.p1.push_back($1);
		i++;
		}
	;

%%