#ifndef __VALUES_H__
#define __VALUES_H__

#ifdef __cplusplus

#include <string>
#include <sstream>

using namespace std;

class Values {
private:
	int numSend;
	int numRecv;
	int sumSrc;
	int sumDest;
	int xorSrc;
	int xorDest;
public:
	Values(int nSend, int nRecv, int sSrc, int sDest, int xSrd, int xDest);
	Values();
	~Values();

	int getNumSend();
	int getNumRecv();
	int getSumSrc();
	int getSumDest();
	int getXorSrc();
	int getXorDest();
	bool match(Values* other);

	string toString();

};

#endif /* __cplusplus */
#endif /*__ITER_H__*/