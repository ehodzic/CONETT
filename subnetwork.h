#ifndef __SUBNETWORK_H_INCLUDED__
#define __SUBNETWORK_H_INCLUDED__
// #include <unordered_map>
#include <map>
#include "bitmask.h"

class SubnetworkEntry {
	// std::unordered_map<int, llu> nodes;	// indexed by node indices, the entry is the alteration bitmask
	std::map<int, llu> nodes;	// indexed by node indices, the entry is the alteration bitmask
	Bitmask * samples;
	bool valid;
	std::string name;

public:
	SubnetworkEntry();
	SubnetworkEntry(int totalNumSamples);
	SubnetworkEntry(const Bitmask * samples);
	SubnetworkEntry(const Bitmask & samples);
	SubnetworkEntry(const SubnetworkEntry & Q);
	// SubnetworkEntry(const std::unordered_map<int, llu> & nodes, const Bitmask & samples);
	SubnetworkEntry(const std::map<int, llu> & nodes, const Bitmask & samples);
	~SubnetworkEntry();

	void operator=(const SubnetworkEntry & Q);

	int getNumSamples() const { return this->samples->getSize(); }
	int getSize() const { return nodes.size(); }
	bool isValid() { return valid; }
	std::string getName() const { return name; }
	// const std::unordered_map<int, llu> & getNodesRef() const { return nodes; }	// Returns a reference to constant nodes map
	const std::map<int, llu> & getNodesRef() const { return nodes; }	// Returns a reference to constant nodes map
	const Bitmask * getSamplesPtr() const { return samples; }	// Returns a pointer to constant samples bitmask

	void addNode(int nodeIdx, int altIdx);
	void setInvalid() { valid = false; }
	void setVaid() { valid = true; }
	void setName(const char * setName) { name = std::string(setName); }
	int removeSamplesWithMoreThanOneAlteration(Bitmask *** geneSampleAlterations);
};

class SubnetworkMask {
	llu * mask;
	int len;

public:
	SubnetworkMask();
	SubnetworkMask(const std::map<int, llu> & nodes, const int bitsInNode, const int bitsInAlt);
	SubnetworkMask(const SubnetworkMask & Q);
	~SubnetworkMask();

	void operator=(const SubnetworkMask & Q);
	bool operator==(const SubnetworkMask & Q) const;
	bool operator<(const SubnetworkMask & Q) const;
	bool operator>(const SubnetworkMask & Q) const;

	std::size_t getHash() const {
		std::size_t res = 17;
		for (int i = 0; i < len; i++) {
			res = res * 31 + std::hash<llu>()(mask[i]);
		}
		return res;
	}
};

struct SubnetworkMaskHasher {
	std::size_t operator()(const SubnetworkMask & Q) const {
		return Q.getHash();
	}
};

#endif