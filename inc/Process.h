#ifndef __PROCESS_H__
#define __PROCESS_H__

#ifdef __cplusplus

#include <queue>

using namespace std;

typedef queue<int> RecvQueue;

typedef struct {
	int noDlks;
	int preRemove;
	RecvQueue priorRecvs;
} Process;



#endif /* __cplusplus */

#endif /* __PROCESS_H__ */