#include "Loop.h"

Loop::Loop():
head(NULL),
tail(NULL)
{}

Loop::~Loop() {}

void Loop::addIter(Iter* iter, int rank) {
	if (!iter) return;
	if (head == NULL) {
		printf("Rank %d : Add First Iter\n", rank);
		head = tail = iter;
		//cout<< "Rank " << rank << ": ";
		//print();
	} else {
		printf("Rank %d : Add Not first Iter\n", rank);
		tail->next = iter;
		iter->prev = tail;
		tail = iter;
	}
}

void Loop::appendIter(Iter* iter, int rank) {
	printf("Rank %d : Enter append Iter\n", rank);
	
	Iter* it = head;
	Values* values = iter->getValues();
	int match_flag = 0;
	
	while(it) {
		printf("Rank %d : Enter while\n", rank);
		//cout<< "Rank " << rank << ": " << it->toString();
		if(it->getValues()->match(values)) {

			match_flag = 1;
			break;
		}
		
		it = it->next;
	}
	printf("Rank %d : Out while\n", rank);

	if(match_flag) {
		printf("Rank %d : Enter match\n", rank);
		it->addIterCount(iter->getIterAt(0));
	
	} else {
		printf("Rank %d : Enter not match\n", rank);
		addIter(iter, rank);
	}
}

void Loop::printLoop(const char* filename, int iter, double time, int rank) {
	ofstream tracefile;
	tracefile.open(filename, ofstream::out | ofstream::app);
	if(!tracefile ){
		cerr << "error: unable to open trace file: " << filename;
		return;
	}
	tracefile << "[" << iter << "]" << "\n";
	Iter* it = head;
	while(it) {
		tracefile << it->toString() << "\n";
		it = it->next;
	}

	if (time != 0){ 
		tracefile  << "!" << time << "\n";
	}
	tracefile.close();
}

void Loop::print() {
	Iter* it = head;
	while(it) {
		cout << it->toString() << "\n";
		it = it->next;
	}
}

void appendIter(Loop* loop, Iter* iter, int rank) {
	loop->appendIter(iter, rank);
}

void printLoop(Loop* loop, const char* filename, int iter, double time, int rank) {
	loop->printLoop(filename, iter, time, rank);
}
