#include "Iter.h"

Iter::Iter(int iter, int nSend, int nRecv, int sSrc, int sDest, int xSrc, int xDest):
prev(NULL),
next(NULL)
{
	iters.push_back(iter);
	values = new Values(nSend, nRecv, sSrc, sDest, xSrc, xDest);
}

Iter::Iter(): 
prev(NULL),
next(NULL) 
{}

Iter::~Iter(){}

Values* Iter::getValues() {
	return values;
}

void Iter::addIterCount(int iter) {
	iters.push_back(iter);
}

int Iter::getIterAt(int i) {
	if (i < 0 || i >= (int)iters.size()) return -1;
	return iters[i];
}

string Iter::toString() {
	ostringstream rtn;
	rtn << "( ";
	for(unsigned i = 0; i < iters.size(); i++) {
		rtn << iters[i] << " ";
	}
	rtn << ")\n";
	
	rtn << values->toString();
	return rtn.str();
}

void createIter(Iter** it, int iter, int nSend, int nRecv, int sSrc, int sDest, int xSrc, int xDest) {
	(*it) = new Iter(iter, nSend, nRecv, sSrc, sDest, xSrc, xDest);
}
