#include "Misc.h"

/* Global Variable */
/*extern int numSend;	//The number of Send Event on each process
extern int numRecv;	//the number of Recv Event on each process
extern int Drank;		//The rank of the current process
extern int sumSrc;		//The xor(src) of all Recv Event on each process
extern int xorSrc;		//The xor(src) of all Recv Event on all processes
extern int sumDest;	//The xor(dest) of all Send Event on each process
extern int xorDest;		//The xor(dest) of all Send Event on all processes 
Loop *loop;
Iter* iter;*/

extern int enabled;

void beginFor() {
	// loop = new Loop();
}
void beginIter() {
	/*numSend = 0;
	numRecv = 0;
	sumSrc = 0;
	sumDest = 0;
	xorSrc = 0;
	xorDest = 0;*/
}

void endIter(int i, int rank) {
	// createIter(&iter, i, numSend, numRecv, sumSrc, sumDest, xorSrc, xorDest);
	// appendIter(loop, iter, rank);
	/*printf("Rank %d : Iter %d : numSend = %d  numRecv = %d sumSrc = %d sumDest = %d xorSrc = %d xorDest = %d\n",rank , i, numSend, numRecv, sumSrc, sumDest, xorSrc, xorDest);*/
}

void endFor(int iter, double time, int rank) {
	/*stringstream temp;
        temp << "traces/" << Drank;
	printLoop(loop, temp.str().c_str(), iter ,time, rank);*/

	/*if (Drank == 0) {
	ofstream iterfile;
	iterfile.open ("traces/iter", ofstream::out | ofstream::app);
	if(!iterfile ){
		cerr << "error: unable to open iteration file ";
		return;
	}
	iterfile << iter << "\n";
	iterfile.close();
	}*/
}

void beginning() {
	enabled = 1;
}

void ending() {
	enabled = 0;
}
