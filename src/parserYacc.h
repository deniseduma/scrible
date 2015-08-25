#define NUMBER 257
#define STRING 258
#define TO 259
#define AVERAGE 260
#define CORREC 261
#define KMERCOUNT 262
typedef union{
int num;
char var[255];
} YYSTYPE;
extern YYSTYPE yylval;
