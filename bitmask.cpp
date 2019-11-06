#include "bitmask.h"

Bitmask::Bitmask(int maxSize) : maxSize(maxSize) {
	size = 0;
	len = ceil(maxSize / 64.0);
	bits = new llu [len];
	memset(bits, 0, sizeof(llu)*len);
}

Bitmask::Bitmask(const Bitmask & Q) {
	maxSize = Q.maxSize;
	size = Q.size;
	len = Q.len;
	bits = new llu [len];
	for (int i = 0; i < len; i++)
		bits[i] = Q.bits[i];
}

Bitmask::Bitmask(const Bitmask * Q) {
	maxSize = Q->maxSize;
	size = Q->size;
	len = Q->len;
	bits = new llu [len];
	for (int i = 0; i < len; i++)
		bits[i] = Q -> bits[i];
}

Bitmask::~Bitmask() {
	if (bits)
		delete bits;
}

bool Bitmask::operator== (const Bitmask & Q) const {
	if (size != Q.size) return false;
	int minLen = (len < Q.len) ? len : Q.len;
	for (int i = 0; i < minLen; i++) if (bits[i] != Q.bits[i]) return false;
	for (int i = minLen; i < len; i++) if (bits[i]) return false;
	for (int i = minLen; i < Q.len; i++) if (Q.bits[i]) return false;
	return true;
}

void Bitmask::operator= (const Bitmask & Q) {
	if (maxSize != Q.maxSize) {
		fprintf(stderr, "< Error in Bitmask::operator= > Bitmasks do not have the same length.\n");
		exit(0);
	}
	maxSize = Q.maxSize;
	size = Q.size;
	len = Q.len;
	for (int i = 0; i < len; i++)
		bits[i] = Q.bits[i];
}

void Bitmask::copylluBitmask(const llu x) {
	len = 1;
	bits[0] = x;
	size = __builtin_popcountll(x);
}

void Bitmask::setBit(int pos, bool val) {
	int idx = pos / 64;
	if (idx >= len) {
		fprintf(stderr, "< Error in Bitmask::setBit > Cannot assign bit because bitmask is too small. pos: %d | len: %d | idx: %d | maxSize: %d\n", pos, len, idx, maxSize);
		exit(0);
	}
	int bitIdx = pos % 64;
	bool oldVal = getBit(pos);
	if (oldVal ^ val) {	// the bit is about to get changed
		bits[idx] ^= llu(1) << bitIdx;
		if (val) size++;
		else size--;
	}
}

bool Bitmask::getBit(int pos) const {
	int idx = pos / 64;
	if (idx >= len) {
		fprintf(stderr, "< Error in Bitmask::getBit > Cannot return bit because bitmask is too small.\n");
		exit(0);
	}
	int bitIdx = pos % 64;
	return bits[idx] & (llu(1) << bitIdx);
}

int Bitmask::getPositionOfFirstSetBit() const {
	if (size == 0) {
		fprintf(stderr, "< Error in Bitmask::getPositionOfFirstSetBit > Cannot return position of first set bit because bitmask is empty.\n");
		exit(0);
	}
	for (int i = 0; i < len; i++) {
		if (bits[i]) return i*64 + __builtin_ctzll(bits[i]);
	}
}

int Bitmask::extractLowestOrderSetBitIndex() {
	if (size == 0) {
		fprintf(stderr, "< Error in Bitmask::extractLowestOrderSetBitIndex > Cannot extract first set bit because bitmask is empty.\n");
		exit(0);
	}
	for (int i = 0; i < len; i++) {
		if (bits[i]) {
			int pos = i*64 + __builtin_ctzll(bits[i]);
			setBit(pos, 0);
			return pos;
		}
	}
}

void Bitmask::negateSelf() {
	size = 0;
	for (int i = 0; i < len; i++) {
		bits[i] = ~bits[i];
		size += __builtin_popcountll(bits[i]);
	}
}

void Bitmask::andWith(const Bitmask & Q) {
	if (maxSize != Q.maxSize) {
		fprintf(stderr, "< Error in Bitmask::andWith > Bitmasks do not have the same length.\n");
		exit(0);
	}
	size = 0;
	for (int i = 0; i < len; i++) {
		bits[i] &= Q.bits[i];
		size += __builtin_popcountll(bits[i]);
	}
}

void Bitmask::andWith(const Bitmask * Q) {
	if (maxSize != Q->maxSize) {
		fprintf(stderr, "< Error in Bitmask::andWith > Bitmasks do not have the same length.\n");
		exit(0);
	}
	size = 0;
	for (int i = 0; i < len; i++) {
		bits[i] &= Q->bits[i];
		size += __builtin_popcountll(bits[i]);
	}
}

void Bitmask::xorWith(const Bitmask & Q) {
	if (maxSize != Q.maxSize) {
		fprintf(stderr, "< Error in Bitmask::xorWith > Bitmasks do not have the same length.\n");
		exit(0);
	}
	size = 0;
	for (int i = 0; i < len; i++) {
		bits[i] ^= Q.bits[i];
		size += __builtin_popcountll(bits[i]);
	}
}

void Bitmask::xorWith(const Bitmask * Q) {
	if (maxSize != Q->maxSize) {
		fprintf(stderr, "< Error in Bitmask::xorWith > Bitmasks do not have the same length.\n");
		exit(0);
	}
	size = 0;
	for (int i = 0; i < len; i++) {
		bits[i] ^= Q->bits[i];
		size += __builtin_popcountll(bits[i]);
	}
}

void Bitmask::orWith(const Bitmask & Q) {
	if (maxSize != Q.maxSize) {
		fprintf(stderr, "< Error in Bitmask::orWith > Bitmasks do not have the same length.\n");
		exit(0);
	}
	size = 0;
	for (int i = 0; i < len; i++) {
		bits[i] |= Q.bits[i];
		size += __builtin_popcountll(bits[i]);
	}
}

void Bitmask::orWith(const Bitmask * Q) {
	if (maxSize != Q->maxSize) {
		fprintf(stderr, "< Error in Bitmask::orWith > Bitmasks do not have the same length.\n");
		exit(0);
	}
	size = 0;
	for (int i = 0; i < len; i++) {
		bits[i] |= Q->bits[i];
		size += __builtin_popcountll(bits[i]);
	}
}

void Bitmask::clear() {
	size = 0;
	for (int i = 0; i < len; i++) {
		bits[i] = 0;
	}
}