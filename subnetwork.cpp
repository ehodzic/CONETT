#include "subnetwork.h"

SubnetworkEntry::SubnetworkEntry() : samples(nullptr), valid(true) {}

SubnetworkEntry::SubnetworkEntry(int totalNumSamples) : valid(true) {
	samples = new Bitmask(totalNumSamples);
}

SubnetworkEntry::SubnetworkEntry(const Bitmask & samples) : valid(true) {
	this->samples = new Bitmask(samples);
}

SubnetworkEntry::SubnetworkEntry(const Bitmask * samples) : valid(true) {
	this->samples = new Bitmask(samples);
}

SubnetworkEntry::SubnetworkEntry(const SubnetworkEntry & Q) : valid(Q.valid) {
	nodes = Q.nodes;
	samples = new Bitmask(Q.samples);
	name = Q.name;
}

// SubnetworkEntry::SubnetworkEntry(const std::unordered_map<int, llu> & nodes, const Bitmask & samples) : valid(true) {
SubnetworkEntry::SubnetworkEntry(const std::map<int, llu> & nodes, const Bitmask & samples) : valid(true) {
	this->nodes = nodes;
	this->samples = new Bitmask(samples);
}

SubnetworkEntry::~SubnetworkEntry(){
	if (samples) {
		delete samples;
		samples = nullptr;
	}
}

void SubnetworkEntry::operator=(const SubnetworkEntry & Q) {
	valid = Q.valid;
	nodes = Q.nodes;
	samples = new Bitmask(Q.samples);
	name = Q.name;
}

void SubnetworkEntry::addNode(int nodeIdx, int altIdx) {
	llu mask = (llu(1) << altIdx);
	if (!nodes.count(nodeIdx))
		nodes[nodeIdx] = mask;
	else 
		nodes[nodeIdx] |= mask;
}

int SubnetworkEntry::removeSamplesWithMoreThanOneAlteration(Bitmask *** geneSampleAlterations) {
	Bitmask tempmask(this->samples);
	int numRemovedSamples = 0;
	while (tempmask.getSize()) {
		int sampleIdx = tempmask.extractLowestOrderSetBitIndex();
		int numAlteredGenes = 0;
		for (auto & it: nodes) {
			int nodeIdx = it.first;
			Bitmask alterationMask(64);
			alterationMask.copylluBitmask(it.second);
			bool geneAltered = false;
			while (alterationMask.getSize()) {
				int altIdx = alterationMask.extractLowestOrderSetBitIndex();
				if (geneSampleAlterations[altIdx][nodeIdx] == nullptr) {
					fprintf(stderr, "< Error in SubnetworkEntry::removeSamplesWithMoreThanOneAlteration > There is a node without a colour in the constructed ME subnetwork.\n");
					exit(0);
				}
				geneAltered |= geneSampleAlterations[altIdx][nodeIdx]->getBit(sampleIdx);
			}
			numAlteredGenes += geneAltered;
		}
		if (numAlteredGenes > 1) {
			samples->setBit(sampleIdx, 0);
			numRemovedSamples++;
		}
	}
	return numRemovedSamples;
}

SubnetworkMask::SubnetworkMask() : len(0), mask(nullptr) {}

SubnetworkMask::SubnetworkMask(const std::map<int, llu> & nodes, const int bitsInNode, const int bitsInAlt) {
	int size = nodes.size();
	int totalBits = bitsInNode + bitsInAlt;
	len = (size * totalBits) / 64 + ((size * totalBits) % 64 != 0);
	mask = new llu [len];
	memset(mask, 0, sizeof(mask[0]) * len);
	int i = 0;
	for (auto const & it : nodes) {
		int nodeIdx = it.first;
		llu altMask = it.second;
		llu nodeMask = (llu(nodeIdx) << bitsInAlt) + altMask;
		int startBit = i * totalBits;
		int startBucket = startBit / 64;
		startBit %= 64;
		int endBit = (i + 1) * totalBits - 1;
		int endBucket = endBit / 64;
		endBit %= 64;

		mask[startBucket] |= nodeMask << startBit;
		if (startBucket < endBucket)
			mask[endBucket] |= nodeMask >> (totalBits - endBit - 1);

		i++;
	}
}

SubnetworkMask::SubnetworkMask(const SubnetworkMask & Q) {
	len = Q.len;
	mask = new llu [len];
	for (int i = 0; i < len; i++) {
		mask[i] = Q.mask[i];
	}
}

void SubnetworkMask::operator=(const SubnetworkMask & Q) {
	if (len != Q.len) {
		delete mask;
		mask = new llu [Q.len];
	}
	len = Q.len;
	for (int i = 0; i < len; i++) {
		mask[i] = Q.mask[i];
	}
}

bool SubnetworkMask::operator==(const SubnetworkMask & Q) const {
	if (len != Q.len)
		return false;
	for (int i = 0; i < len; i++) {
		if (mask[i] != Q.mask[i])
			return false;
	}
	return true;
}

bool SubnetworkMask::operator<(const SubnetworkMask & Q) const {
	if (len != Q.len)
		return (len < Q.len);
	for (int i = 0; i < len; i++) {
		if (mask[i] != Q.mask[i])
			return (mask[i] < Q.mask[i]);
	}
	return false;
}

bool SubnetworkMask::operator>(const SubnetworkMask & Q) const {
	if (len != Q.len)
		return (len > Q.len);
	for (int i = 0; i < len; i++) {
		if (mask[i] != Q.mask[i])
			return (mask[i] > Q.mask[i]);
	}
	return false;
}

SubnetworkMask::~SubnetworkMask(){
	if (mask != nullptr) {
		delete mask;
		mask = nullptr;
	}
}