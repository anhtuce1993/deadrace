#ifndef __LOOP_H__
#define __LOOP_H__



#ifdef __cplusplus

#include <fstream>
#include <iostream>

#include "Iter.h"

using namespace std;

class Loop {
private:
	Iter* head;
	Iter* tail;
public:
	Loop();
	~Loop();

	void appendIter(Iter* iter, int rank);

	void printLoop(const char* filename, int iter, double time, int rank);

	void addIter(Iter* iter, int rank);

	void print();

};

extern "C" {
#endif /* __cplusplus */


typedef struct Loop Loop;

void appendIter(Loop* loop, Iter* iter, int rank);

void printLoop(Loop* loop, const char* filename, int iter, double time, int rank);

#ifdef __cplusplus 
}
#endif /* __cplusplus */

#endif /*__LOOP_H__*/
