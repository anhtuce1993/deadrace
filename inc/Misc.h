#ifndef __MISC_H__
#define __MISC_H__

#include "Iter.h"
#include "Loop.h"


/* Global Variable */
/*extern int numSend;	//The number of Send Event on each process
extern int numRecv;	//the number of Recv Event on each process
extern int Drank;		//The rank of the current process
extern int sumSrc;		//The xor(src) of all Recv Event on each process
extern int xorSrc;		//The xor(src) of all Recv Event on all processes
extern int sumDest;	//The xor(dest) of all Send Event on each process
extern int xorDest;		//The xor(dest) of all Send Event on all processes 
extern Loop *loop;
extern Iter* iter;*/

extern int enabled;

void beginFor();
void beginIter();
void endIter(int i, int rank);
void endFor(int iter, double time, int rank);

extern "C" void beginning();
extern "C" void ending();

/*void beginning();
void ending();*/



#endif /* __MISC_H__ */
