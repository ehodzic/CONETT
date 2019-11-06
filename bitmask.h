#ifndef __BITMASK_H_INCLUDED__
#define __BITMASK_H_INCLUDED__
#include <cstring>
#include <cstdio>
#include <iostream>
#include <cmath>
typedef unsigned long long llu;

struct Bitmask {
	int maxSize;
	int size;
	int len;
	llu * bits;

	Bitmask(int maxSize);
	Bitmask(const Bitmask & Q);
	Bitmask(const Bitmask * Q);
	~Bitmask();

	bool operator== (const Bitmask & Q) const;
	void operator= (const Bitmask & Q);

	bool getBit(int pos) const;
	int getSize() const { return size; }
	int getPositionOfFirstSetBit() const;

	void setBit(int pos, bool val);

	void copylluBitmask(const llu x);
	int extractLowestOrderSetBitIndex();
	void negateSelf();
	void andWith(const Bitmask & Q);
	void andWith(const Bitmask * Q);
	void xorWith(const Bitmask & Q);
	void xorWith(const Bitmask * Q);
	void orWith(const Bitmask & Q);
	void orWith(const Bitmask * Q);
	void clear();
};

struct BitmaskHasher {
	std::size_t operator()(const Bitmask & Q) const {
		using std::size_t;
		using std::hash;
		size_t res = 17;
		for (int i = 0; i < Q.len; i++) {
			res = res * 31 + hash< unsigned long long >()(Q.bits[i]);
		}
		return res;
	}
};

#endif