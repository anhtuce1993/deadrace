#include "Controller.h"

Controller::Controller() {
	increment = 0;
}

Controller::~Controller() {}

/*void Controller::insertRecv(int rank, int clock) {
	// Process pqueue;
	RecvQueue pqueue;
	map<int, RecvQueue >::iterator it = commProcs.find(rank);
	if (it == commProcs.end()) {
		// pqueue.priorRecvs.push(clock);
		// commProcs.insert(pair<int, Process>(rank, pqueue));
		pqueue.push(clock);
		commProcs.insert(pair<int, RecvQueue >(rank, pqueue));
	} else {
		pqueue = it->second;
		// pqueue.priorRecvs.push(clock);
		pqueue.push(clock);
	}
}*/

void Controller::addRootRecv(int lclk, int from) {
	rootRecvs.push_back(from);
}

void Controller::printRootRecvs() {
	printf("\nRoot Receiving List Remains: ");
	for(unsigned i = 0; i < rootRecvs.size(); i++) {
		printf("\n%i ", i+1+increment);
		(rootRecvs[i] == -1) ? printf("MPI_ANY_SOURCE") : printf("%i", rootRecvs[i]);
	}
}

void Controller::printProcRecvs(int src) {
	map<int, Process>::iterator it = commProcs.find(src);
	// map<int, RecvQueue >::iterator it = commProcs.find(src);
	if (it == commProcs.end()) {
		printf("\nReceiving queue of process %i has not been initialized", src);
	} else {
		// Process pqueue = it->second;
		RecvQueue pqueue = it->second.priorRecvs;
		// queue<int> tmpqueue = it->second.priorRecvs;
		// int tmp;
		printf("\nProcess %i : ", src);
		printf("[front = %i back = %i] ", pqueue.front(), pqueue.back());
		// while(!tmpqueue.empty()) {
		// 	tmp = tmpqueue.front();
		// 	printf("%i ", tmp);
		// 	tmpqueue.pop();
		// 	// tmpqueue.push(tmp);
		// }
	}
}

int Controller::minPreRemove(int nProc, int rootProc) {
	/* 0 : there is a process not initialized 
	   > 0 : min receiving iterator */
	int min = -1;
	map<int, Process >::iterator it;
	for (int i = 0; i < nProc; i++) {
		if (i != rootProc) {
			/*printf(" i = %i",i);*/
			it = commProcs.find(i);
			if (it == commProcs.end()) return 0; 
			if (min < 0) {
				min = it->second.preRemove;
				/*printf("min < 0");*/
			}
			else {
				min = (it->second.preRemove < min) ? it->second.preRemove : min; 
				/*printf("min > 0");*/
			}
		}
	}
	return min;
}

void Controller::manipulateCommProcs(int src, int recvlclk, int lclk, FILE* file) {
	// Process pqueue;
	RecvQueue recvs;
	unsigned i;
	// map<int, RecvQueue >::iterator it = commProcs.find(src);
	map<int, Process >::iterator it = commProcs.find(src);
	if (it == commProcs.end()) {
		// receiving queue of src process have not been initialized

		// add receiving queue of src process into map
		// int first = 0;
		/*printf("not initalized");*/
		Process pqueue;
		for (i = recvlclk; i < lclk; i++) {
			if (rootRecvs[i-increment] == -1 || rootRecvs[i-increment] == src) {
				/*if (first == 0) {
					first = 1;	
				} else {
					pqueue.priorRecvs.push(i+1);
				}*/
				// pqueue.priorRecvs.push(i+1);
				recvs.push(i+1);
			}
		}
		/*printf("[front = %i back = %i] ", recvs.front(), recvs.back());*/
		// pqueue.priorRecvs.pop();
		pqueue.preRemove = recvs.front();
		recvs.pop();
		/*printf("[front = %i back = %i] ", recvs.front(), recvs.back());*/
		if (recvs.empty()) printf(" Empty!!! ");
		
		pqueue.priorRecvs = recvs;
		pqueue.noDlks = 0;
		// pqueue.preRemove = 0;
		commProcs.insert(pair<int, Process>(src, pqueue));
		// commProcs.insert(pair<int, RecvQueue >(src, recvs));
	} else {
		/*printf("initialized");*/
		// receiving queue of src process have been initialized
		// pqueue = it->second;
		// queue<int> rqueue = pqueue.priorRecvs;
		recvs = (it->second).priorRecvs;
		int dlks = (it->second).noDlks;
		int preRm = (it->second).preRemove;
		// int max = (recvlclk > rqueue.back()) ? recvlclk + 1 : rqueue.back() + 1;
		// printf("[back = %i] ", recvs.back());
		printf(" preRm : %i ",preRm);
		if (recvs.empty()) {

			/*printf(" empty ");*/
			
			/* In this step, pushing in queue at least one element */
			int t = (preRm > recvlclk) ? preRm : recvlclk;
			for (i = t - increment; i < rootRecvs.size(); i++) {
				if (rootRecvs[i] == -1 || rootRecvs[i] == src)
				// add receiving queue of src process into map
					recvs.push(i+1+increment);
			}
			/* So there is no case of popping being taken on empty queue */
			/*printf("[front = %i back = %i] ", recvs.front(), recvs.back());*/
			it->second.preRemove = recvs.front();
			recvs.pop();
			/*printf("[front = %i back = %i] ", recvs.front(), recvs.back());*/
			if (recvs.empty()) printf(" Empty!!! ");
		} else {
			printf("not empty");
			/* Queue is not empty so it is OK to get back at this step */
			int back = recvs.back(); 
			/*printf(" %i ", back);*/

			/* CASE iterator of recvs remaining in queue is less than receving local clock */
			if (recvlclk > back) {
				/*printf(" case1 ");*/
				while(!recvs.empty()) {
					if (rootRecvs[recvs.front()-1-increment] == src) {
						fprintf(file, "Deadlock happens at RC = %i , RECV( %i )\n", recvs.front(), src);
						printf(" Deadlock ");
						dlks++;
					}
					recvs.pop();
				}
				for (i = back; i < recvlclk; i++) {
					if (rootRecvs[i-increment] == src) {
						fprintf(file, "Deadlock happens at RC = %i , RECV( %i )\n", i + 1, src);
						printf(" Deadlock ");
						dlks++;
					}
				} 
				for (i = recvlclk - increment; i < rootRecvs.size(); i++) {
					if (rootRecvs[i] == -1 || rootRecvs[i] == src)
					/* add receiving queue of src process into map */
						recvs.push(i+1+increment);
				}
				/*printf("[front = %i back = %i] ", recvs.front(), recvs.back());*/
				it->second.preRemove = recvs.front();
				recvs.pop();
				/*printf("[front = %i back = %i] ", recvs.front(), recvs.back());*/
				if (recvs.empty()) printf(" Empty!!! ");
			}
			/* CASE iterator of recvs remaining in queue is greater than receving local clock */ 
			else {
				/*printf(" case2 recvlclk=%i back=%i ", recvlclk, back);*/
				int tmp;
				while (!recvs.empty() && (tmp = recvs.front()) <= recvlclk) {
					/*printf("[front = %i back = %i] ", tmp, recvs.back());*/
					if (rootRecvs[tmp - 1 - increment] == src) {
						fprintf(file, "Deadlock happens at RC = %i , RECV( %i )\n", tmp, src);
						printf(" Deadlock ");
						dlks++;
					}
					recvs.pop();
				}
				/*printf(" EndRemove ");*/
				/*printf("[front = %i back = %i] ", recvs.front(), recvs.back());*/
				for (i = back - increment; i < rootRecvs.size(); i++) {
					if (rootRecvs[i] == -1 || rootRecvs[i] == src) {
					/*add receiving queue of src process into map*/
						// printf("Enter %i", i);
						recvs.push(i+1+increment);
						/*printf("[front = %i back = %i] ", recvs.front(), recvs.back());*/
					}
				}
				/*printf("[front = %i back = %i] ", recvs.front(), recvs.back());*/
				it->second.preRemove = recvs.front();
				recvs.pop();
				/*printf("[front = %i back = %i] ", recvs.front(), recvs.back());*/
				if (recvs.empty()) printf(" Empty!!! ");
			}
		}
		it->second.priorRecvs = recvs;
		it->second.noDlks = dlks;
	}
}

void Controller::checkRemainQueue(int rank, FILE* file) {
	map<int, Process >::iterator it = commProcs.find(rank);
	// map<int, RecvQueue >::iterator it = commProcs.find(rank);
	if (it == commProcs.end()) {
		printf("\nProcess %d check queue: not initalized \n", rank);
		return;
	} else {
		// at MPI_Finalize queue is empty which also mean that all possible recvs is added or just remain non-added recvs (*)
		// so definitely not happen deadlock in this queue causing by recvs standing outside queue
		printf("\nProcess %d check queue: initalized \n", rank);
		RecvQueue recvs = it->second.priorRecvs;
		int dlks = (it->second).noDlks;
		/*fprintf(file, "noDlks = %i", dlks);*/
		// printf("[front = %i back = %i] ", recvs.front(), recvs.back());
		if (!recvs.empty()) {
			int tmp;
			while (!recvs.empty()) {
				tmp = recvs.front();
				if (rootRecvs[tmp - 1 - increment] == rank) {
					fprintf(file, "Deadlock happens at RC = %i , RECV( %i )\n", tmp, rank);
					dlks++;
					// return;
				}
				recvs.pop();
			}
		} /*else {
			// at MPI_Finalize queue is empty which also mean that all possible recvs is added or just remain non-added recvs (*)
			// so definitely not happen deadlock in this queue causing by recvs standing outside queue
		}
		int i;
		int back = recvs.back();
		printf("\nBack of queue : %d", back);
		for (i = back; i < rootRecvs.size(); i++) {
			if (rootRecvs[i] == rank) {
				printf(" Deadlock!!! ");
				return 1;
			}
		}*/
		if (dlks > 0) {
			printf("\tDealock happens. No of deadlocks happening : %i\n", dlks);
		} else {
			printf("\tNo dealock !!!\n");
		}
		return ;
	}
}

void Controller::removeRootRecvs(int toIter) {
	if (toIter > increment) {
		rootRecvs.erase(rootRecvs.begin(), rootRecvs.begin() + toIter - increment);
		increment = toIter;
	}
}

void initController(Controller **controller) {
	(*controller) = new Controller();
}

void addRootRecv(Controller* controller, int lclk, int from) {
	controller->addRootRecv(lclk,from);
}

void printRootRecvs(Controller *controller) {
	controller->printRootRecvs();
}

void manipulateCommProcs(Controller* controller, int src, int recvlclk, int lclk, FILE* file) {
	controller->manipulateCommProcs(src, recvlclk, lclk, file);
}

void printProcRecvs(Controller* controller, int src) {
	controller->printProcRecvs(src);
}

void checkRemainQueue(Controller* controller, int rank, FILE* file) {
	controller->checkRemainQueue(rank,file);
	/*if (controller->checkRemainQueue(rank) == 1) {
		printf("\nDeadlock happens from sending of process %d ", rank);
	} else {
		printf("\nDeadlock does not happens from sending of process %d ", rank);
	}*/
}

int minPreRemove(Controller* controller, int nProc, int rootProc) {
	return controller->minPreRemove(nProc,rootProc);
}

void removeRootRecvs(Controller* controller, int toIter) {
	controller->removeRootRecvs(toIter);
}