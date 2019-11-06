#ifndef __GRAPH_H_INCLUDED__
#define __GRAPH_H_INCLUDED__
#include <unordered_map>
#include <unordered_set>
#include <cstdio>
#include "bitmask.h"

struct Graph {
	int V, E;
	int * NSize;
	int ** N;
	std::string * nodeNames;
	int * chrArm;
	std::unordered_map<std::string, int> nodeIndices;

	Graph();
	~Graph();

	void readDirectedEdges(const char * filename, bool writeOutput = false);
};

struct stringPairHash {
	size_t operator() (const std::pair<std::string, std::string> & Q) const {
		size_t h1 = std::hash<std::string>() (Q.first);
		size_t h2 = std::hash<std::string>() (Q.second);
		return h1 ^ (h2 << 1);
	}
};

#endif