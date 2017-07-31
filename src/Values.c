#include "Values.h"

Values::Values(int nSend, int nRecv, int sSrc, int sDest, int xSrc, int xDest):
numSend(nSend), 
numRecv(nRecv),
sumSrc(sSrc),
sumDest(sDest),
xorSrc(xSrc),
xorDest(xDest) {

}

Values::Values(){

}

Values::~Values(){

}

int Values::getNumSend() {
	return numSend;
}

int Values::getNumRecv() {
	return numRecv;
}

int Values::getSumSrc() {
	return sumSrc;
}

int Values::getSumDest() {
	return sumDest;
}

int Values::getXorSrc() {
	return xorSrc;
}

int Values::getXorDest() {
	return xorDest;
}

bool Values::match(Values *other) {
	return ((numSend == other->getNumSend()) && (numRecv == other->getNumRecv()) && (sumSrc == other->getSumSrc()) && (sumDest == other->getSumDest()
		) && (xorSrc == other->getXorSrc()
		) && (xorDest == other->getXorDest()
		));

	/*return ((numSend - other->getNumSend()) & (numRecv - other->getNumRecv()) & (sumSrc - other->getSumSrc()) & (sumDest - other->getSumDest
		) & (xorSrc - other->getXorSrc()
		) & (xorDest - other->getXorDest()
		) == 0);*/
}

string Values::toString() {
	ostringstream rtn;

	rtn << numSend << "\n";
	rtn << numRecv << "\n";
	rtn << sumSrc << "\n";
	rtn << sumDest << "\n";
	rtn << xorSrc << "\n";
	rtn << xorDest << "\n";

	return rtn.str();
}