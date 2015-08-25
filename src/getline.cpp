#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstring>

using namespace std;

int main() {
char line[5];
char line2[5];
ifstream in;
in.open("getline.txt", ios_base::in);
while (in.good()) {
	line[0]='\0';
	in.getline(line, 5);
	printf("%s %lu, %lu\n", line, strlen(line), sizeof(line));
	//strncpy(line2, line, 4);
	//if (in.gcount()>0)
	//	printf("%s, num chars read %lu\n", line2, strlen(line2));
}
in.close();
return 0;
}
