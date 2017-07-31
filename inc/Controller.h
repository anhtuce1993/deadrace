#ifndef __CONTROLLER_H__
#define __CONTROLLER_H__ 


#include <stdio.h>
#include <stdlib.h>


#ifdef __cplusplus

#include <vector>
#include <map>
#include <queue>

#include "Process.h"


using namespace std;

class Controller {
private:
	// map<int, RecvQueue> commProcs;
	map<int, Process> commProcs;
	vector<int> rootRecvs;
	int increment;
public:
	Controller();
	~Controller();

	/*void insertRecv(int process, int clock);*/

	void addRootRecv(int lclk, int from);

	void printRootRecvs();
	
	void manipulateCommProcs(int src, int recvlclk, int lclk, FILE* file);

	void printProcRecvs(int src);

	void checkRemainQueue(int rank, FILE* file);

	int minPreRemove(int nProc, int rootProc);

	void removeRootRecvs(int toIter);

};

void initController(Controller** controller);

void addRootRecv(Controller* controller, int lclk, int from);

void printRootRecvs(Controller* controller);

void manipulateCommProcs(Controller* controller, int src, int recvlclk, int lclk, FILE* file);

void printProcRecvs(Controller* controller, int src);

void checkRemainQueue(Controller* controller, int rank, FILE *file);

int minPreRemove(Controller* controller, int nProc, int rootProc);

void removeRootRecvs(Controller* controller, int toIter);

#endif /* __cplusplus */

#endif /* __CONTROLLER_H__ */
