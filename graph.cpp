#include "graph.h"

Graph::Graph() : V(0), E(0), NSize(nullptr), N(nullptr), nodeNames(nullptr), chrArm(nullptr) {}

Graph::~Graph() {
	for (int i = 0; i < V; i++) {
		if (N[i] != nullptr) {
			delete N[i];
		}
	}
	delete NSize;
	delete N;
	delete nodeNames;
	delete chrArm;
}

/*
** Reads the input -n parameter as a collection of directed edges - pairs of node names that are separated by whitespace.
** The first node is the source and the second is the destination.
*/
void Graph::readDirectedEdges(const char * filename, bool writeOutput) {
	if (writeOutput)
		fprintf(stderr, "Reading the network... ");
	int timerStart = clock();
	std::unordered_set<std::pair<std::string, std::string>, stringPairHash> uniqueEdges;
	char u[1000], v[1000];
	FILE * fin = NULL;
	if (!(fin = fopen(filename, "r"))) {
		fprintf(stderr, "\n< Error > Cannot open file '%s'. Please make sure the file exists.\n", filename);
		exit(0);
	}
	V = 0;
	while (fscanf(fin, "%s%s", u, v) == 2) {
		char * a = u, * b = v;
		if (strcmp(a, b) == 0) continue;
		std::pair<std::string, std::string> e = std::make_pair(std::string(a), std::string(b));
		if (!uniqueEdges.count(e)) {	// Previously NOT seen directed edge
			uniqueEdges.insert(e);
			if (!nodeIndices.count(e.first)) {
				nodeIndices[e.first] = V;
				V++;
			}
			if (!nodeIndices.count(e.second)) {
				nodeIndices[e.second] = V;
				V++;
			}
		}
	}
	fclose(fin);
	E = uniqueEdges.size();

	NSize = new int [V];
	N = new int * [V];
	nodeNames = new std::string[V];
	chrArm = new int [V];
	memset(NSize, 0, sizeof(NSize[0]) * V);
	memset(chrArm, 0, sizeof(chrArm[0]) * V);

	for (auto e: uniqueEdges) {
		int idx1 = nodeIndices[e.first];
		if (nodeNames[idx1].size() == 0) nodeNames[idx1] = e.first;

		int idx2 = nodeIndices[e.second];
		if (nodeNames[idx2].size() == 0) nodeNames[idx2] = e.second;

		NSize[idx1]++;
	}

	for (int i = 0; i < V; i++) {
		N[i] = new int[ NSize[i] ];
	}

	int * tempNSize = new int [V];
	memset(tempNSize, 0, sizeof(tempNSize[0]) * V );

	for (auto e: uniqueEdges) {
		int idx1 = nodeIndices[e.first];
		int idx2 = nodeIndices[e.second];
		N[ idx1 ][ tempNSize[idx1] ] = idx2;
		tempNSize[idx1]++;
	}
	delete [] tempNSize;
	if (writeOutput)
		fprintf(stderr, "done. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC);
	if (writeOutput)
		fprintf(stderr, "\tInput network contains %d nodes and %d directed edges (self-loops are removed).\n", V, E/2);
}