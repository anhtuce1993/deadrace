#ifndef __ITER_H__
#define __ITER_H__

#ifdef __cplusplus

#include <vector>
#include <string>
#include <sstream>

#include "Values.h" 

using namespace std;

class Iter {
private:
	vector<int> iters;
	Values* values; 
public:
	Iter* prev;
	Iter* next;

	Iter(int iter, int nSend, int nRecv, int sSrc, int sDest, int xSrc, int xDest);
	Iter();
	~Iter();

	void addIterCount(int iter);

	int getIterAt(int i);

	Values* getValues();

	string toString();

};

extern "C" {
#endif /* __cplusplus */

typedef struct Iter Iter;

void createIter(Iter** it, int iter, int nSend, int nRecv, int sSrc, int sDest, int xSrc, int xDest);

#ifdef __cplusplus 
}
#endif /* __cplusplus */


#endif /* __ITER_H__ */