//OPT: -lpthread
#include<thread>
#include<sstream>
#include<iostream>
#include<vector>

using namespace std;

void print100(int threadnum) {
	for(int i=0;i<100;i++)
		cout << "thread " << threadnum << ", line " << i << endl;
	}

void print100plus(int threadnum, int &x) {
	for(int i=0;i<100;i++)
		cout << "thread " << threadnum << " (" << this_thread::get_id() << ") " << ", line " << i << endl;
	x++;
}

int main(int argc, char **argv) {
thread t1([]{ for(int i=0;i<1000000;i++) cout << "thread 1, line " << i << endl; });
thread t2([]{ for(int i=0;i<1000000;i++) cout << "thread 2, line " << i << endl; });

//thread t1(print100,1);
//thread t2(print100,2);
												cout << "hardware reported concurrency support = " << thread::hardware_concurrency() << endl;

//int xx;
//thread t1(print100plus,1,std::ref(xx));
//thread t2(print100plus,2,std::ref(xx));

t1.join();
t2.join();

//cout << "xx = " << xx << endl;
}
