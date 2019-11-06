#include <cstdio>
#include <cstring>
#include <ctime>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <random>
#include "bitmask.h"
#include "subnetwork.h"
#include "graph.h"
#include "entry.h"
#include "heap.h"
#include "gurobi_c++.h"
using namespace std;

typedef unsigned long long llu;
typedef long long ll;

// Fancy Output Functions
void printHeader(const char * text) {
	int n = strlen(text);
	fprintf( stderr, "\n" );
	for (int i = 0; i < n + 4; i++) fprintf(stderr, "*");
	fprintf(stderr, "\n* %s *\n", text);
	for (int i = 0; i < n + 4; i++) fprintf(stderr, "*");
	fprintf(stderr, "\n");
}

bool cmpPair(const pair<int, int> & A, const pair<int, int> & B) {
	// Compare based on second number
	return A.second < B.second;
}

void readAlterationPhylogeny(const char * filename, Entry & patients, Entry & nodes, Entry & alterations, Bitmask *** & geneSampleAlterations, Graph & G, unordered_map<string, Bitmask> & edgePatients, const char * outputFolder) {
	fprintf(stderr, "%s\n", filename);
	fprintf(stderr, "Reading the alteration phylogeny... ");
	int timerStart = clock();
	FILE * fin = NULL;
	if (!(fin = fopen(filename, "r"))) {
		fprintf(stderr, "\n< Error in readAlterationPhylogeny > Cannot open file '%s'. Please make sure the file exists.\n", filename);
		exit(0);
	}
	char patient[1000];
	char node1[1000];
	char alt1[1000];
	char node2[1000];
	char alt2[1000];
	unordered_set<pair<string, string>, stringPairHash> uniqueEdges;
	// fscanf(fin, "%*s%*s%*s\n");	// Skip the first row (header)
	while (fscanf(fin, "%s\t%s\t%s\t%s\t%s", patient, node1, alt1, node2, alt2) == 5) {
		if (!patients.indices.count(string(patient))) {
			int idx = patients.getSize();
			patients.indices[string(patient)] = idx;
		}
		if (!nodes.indices.count(string(node1))) {
			int idx = nodes.getSize();
			nodes.indices[string(node1)] = idx;
		}
		if (!nodes.indices.count(string(node2))) {
			int idx = nodes.getSize();
			nodes.indices[string(node2)] = idx;
		}
		// if (string(alt1) != "-" && !alterations.indices.count(string(alt1))) {
		if (!alterations.indices.count(string(alt1))) {
			int idx = alterations.getSize();
			alterations.indices[string(alt1)] = idx;
		}
		// if (string(alt2) != "-" && !alterations.indices.count(string(alt2))) {
		if (!alterations.indices.count(string(alt2))) {
			int idx = alterations.getSize();
			alterations.indices[string(alt2)] = idx;
		}
		pair<string, string> e = make_pair(string(node1), string(node2));
		uniqueEdges.insert(e);
	}
	patients.names = new string[patients.getSize()];
	nodes.names = new string[nodes.getSize()];
	alterations.names = new string[alterations.getSize()];
	for (auto it : patients.indices) {
		patients.names[it.second] = it.first;
	}
	for (auto it : nodes.indices) {
		nodes.names[it.second] = it.first;
	}
	for (auto it : alterations.indices) {
		alterations.names[it.second] = it.first;
	}
	G.E = uniqueEdges.size();
	G.V = nodes.indices.size();
	G.NSize = new int [G.V];
	G.N = new int * [G.V];
	memset(G.NSize, 0, sizeof(G.NSize[0]) * G.V);
	for (auto e: uniqueEdges) {
		int idx1 = nodes.indices[e.first];
		G.NSize[idx1]++;
	}
	for (int i = 0; i < G.V; i++) {
		G.N[i] = new int[ G.NSize[i] ];
	}
	int * tempNSize = new int [G.V];
	memset(tempNSize, 0, sizeof(tempNSize[0]) * G.V );
	for (auto e: uniqueEdges) {
		int idx1 = nodes.indices[e.first];
		int idx2 = nodes.indices[e.second];
		G.N[ idx1 ][ tempNSize[idx1] ] = idx2;
		tempNSize[idx1]++;
	}
	delete [] tempNSize;
	rewind(fin);

	// geneAlterationBitmasks = new unordered_map<int, llu> [ G.V ];
	geneSampleAlterations = new Bitmask ** [alterations.getSize()];
	for (int i = 0; i < alterations.getSize(); i++) {	// For every alteration type
		geneSampleAlterations[i] = new Bitmask * [G.V];
		for (int j = 0; j < nodes.getSize(); j++)	// For every gene
			geneSampleAlterations[i][j] = nullptr;	// We have a binary sample alteration vector
	}
	while (fscanf(fin, "%s\t%s\t%s\t%s\t%s", patient, node1, alt1, node2, alt2) == 5) {
		int patientIndex = patients.indices[patient];
		int node1Index = nodes.indices[string(node1)];
		int node2Index = nodes.indices[string(node2)];
		// if ((string(alt1) == "CNLOSS" || string(alt1) == "CNGAIN") && (string(alt2) == "CNLOSS" || string(alt2) == "CNGAIN")) {
		// 	string label1 = string(node1).substr(5);
		// 	string label2 = string(node2).substr(5);
		// 	int num1, num2;
		// 	sscanf(label1.c_str(), "%d", &num1);
		// 	sscanf(label2.c_str(), "%d", &num2);
		// 	if (num1 == num2) {
		// 		fprintf(stderr, "%s\t%s\t%s\t%s\t%s\n", patient, node1, alt1, node2, alt2);
		// 	}
		// }
		if (alterations.indices.count(string(alt1))) {
			int alteration1Index = alterations.indices[alt1];
			// if ( !geneAlterationBitmasks[node1Index].count(patientIndex) )
			// 	geneAlterationBitmasks[node1Index][patientIndex] = 0;
			// llu bitMask1 = (llu(1) << alteration1Index);
			// bitMask1 |= geneAlterationBitmasks[node1Index][patientIndex];
			// geneAlterationBitmasks[node1Index][patientIndex] = bitMask1;
			if (geneSampleAlterations[alteration1Index][node1Index] == nullptr)
				geneSampleAlterations[alteration1Index][node1Index] = new Bitmask(patients.getSize());
			geneSampleAlterations[alteration1Index][node1Index]->setBit(patientIndex, 1);
		}
		if (alterations.indices.count(string(alt2))) {
			int alteration2Index = alterations.indices[alt2];
			// if ( !geneAlterationBitmasks[node2Index].count(patientIndex) )
			// 	geneAlterationBitmasks[node2Index][patientIndex] = 0;
			// llu bitMask2 = (llu(1) << alteration2Index);
			// bitMask2 |= geneAlterationBitmasks[node2Index][patientIndex];
			// geneAlterationBitmasks[node2Index][patientIndex] = bitMask2;
			if (geneSampleAlterations[alteration2Index][node2Index] == nullptr)
				geneSampleAlterations[alteration2Index][node2Index] = new Bitmask(patients.getSize());
			geneSampleAlterations[alteration2Index][node2Index]->setBit(patientIndex, 1);
		}
		string edgeName = string(node1) + " " + string(alt1) + " " + string(node2) + " " + string(alt2);
		edgePatients.emplace(edgeName, patients.getSize());
		edgePatients.find(edgeName)->second.setBit(patientIndex, 1);
		// if (edgePatients.find(edgeName)->second.getSize() > 1)
		// 	fprintf(stderr, "%s\t%d\n", edgeName.c_str(), edgePatients.find(edgeName)->second.getSize());
	};

	fclose(fin);
	fprintf(stderr, "done. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC);
	fprintf(stderr, "\tThere are %d patients, with a total of %d nodes (altered genes or CNA events), harboring %d different alterations.\n", patients.getSize(), nodes.getSize(), alterations.getSize());
	fprintf(stderr, "Maximum recurrence for each alteration:\n");
	for (int altIdx = 0; altIdx < alterations.getSize(); altIdx++) {
		FILE * foutAlt = fopen((string(outputFolder) + "/distribution_" + alterations.names[altIdx] + ".txt").c_str(), "w");
		fprintf(foutAlt, "SampleRecurrence\tNodeDegree\n");
		int maxIdx = -1;
		for (int nodeIdx = 0; nodeIdx < nodes.getSize(); nodeIdx++) {
			if (geneSampleAlterations[altIdx][nodeIdx] != nullptr) {
				// fprintf(foutAlt, "%d\t%d\n", geneSampleAlterations[altIdx][nodeIdx]->getSize(), G.NSize[nodeIdx]);
				fprintf(foutAlt, "%s\t%d\t%d\n", nodes.names[nodeIdx].c_str(), geneSampleAlterations[altIdx][nodeIdx]->getSize(), G.NSize[nodeIdx]);
				if (maxIdx == -1 || geneSampleAlterations[altIdx][nodeIdx]->getSize() > geneSampleAlterations[altIdx][maxIdx]->getSize())
					maxIdx = nodeIdx;
			}
		}
		fprintf(stderr, "\t%s\t%d\n", alterations.names[altIdx].c_str(), geneSampleAlterations[altIdx][maxIdx]->getSize());
		fclose(foutAlt);
	}
}

void constructSingleNodeSubnetworks(int minDepth, SubnetworkEntry * & output_CandidateSubnetworks, llu & output_NumSubnetworks, Entry & patients, Entry & nodes, Entry & alterations, Bitmask *** geneSampleAlterations, bool writeOutputFlag = true) {
	int timerStart;
	vector<SubnetworkEntry> candidateSubnetworks;
	
	// Construct single-node subnetworks
	timerStart = clock();
	// Bitmask P(patients.getSize());
	for (int g1 = 0, lastProg = 0; g1 < nodes.getSize(); g1++) {
		if (writeOutputFlag) {
			int progress = 1000 * double(g1 + 1) / double(nodes.getSize());
			if (progress > lastProg) {
				// fprintf(stderr, "\r%.1lf%%", double(progress)/10.0);
				lastProg = progress;
			}
		}
		if (nodes.names[g1] == "GL") continue;	// Skipping GL seed
		// fprintf(stderr, "%s\t", nodes.names[g1].c_str());
		for (int i1 = 0; i1 < alterations.getSize(); i1++) {
			// if (alterations.names[i1] != "SNV") continue;	// Skipping all seeds except SNVs
			// if (alterations.names[i1] != "DEL") continue;	// Skipping all seeds except deletions
			// if (alterations.names[i1] != "SDNV" && alterations.names[i1] != "INDEL") continue;	// Skipping all seeds except S/DNVs and INDELs
			// if (alterations.names[i1] != "SNV" && alterations.names[i1] != "DNV" && alterations.names[i1] != "INS" && alterations.names[i1] != "DEL") continue;	// Skipping all seeds except S/DNVs and INDELs
			// if (alterations.names[i1].substr(0, 4) == "EXPR") continue;	// Skipping expression outlier seeds
			// if (nodes.names[g1] != "CDKN2A" || alterations.names[i1] != "DEL") continue;	// Skipping expression outlier seeds
			// if (g1 == nodes.indices["ACTR5"] && geneSampleAlterations[i1][g1] != nullptr)
			// 		fprintf(stderr, "%d\n", geneSampleAlterations[i1][g1]->getSize());
			if (geneSampleAlterations[i1][g1] != nullptr && geneSampleAlterations[i1][g1]->getSize() >= minDepth) {	// If g1 has alterations of type i1 and of sufficient depth
				// fprintf(stderr, "%d\t", geneSampleAlterations[i1][g1]->getSize());
				SubnetworkEntry newEntry(geneSampleAlterations[i1][g1]);
				newEntry.addNode(g1, i1);
				candidateSubnetworks.push_back(newEntry);
			}
		}
	}
	// fprintf(stderr, "\r\tDone. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC);
	
	output_CandidateSubnetworks = new SubnetworkEntry[candidateSubnetworks.size()];
	output_NumSubnetworks = 0;
	for (int i = 0; i < candidateSubnetworks.size(); i++) {
		output_CandidateSubnetworks[output_NumSubnetworks] = candidateSubnetworks[i];
		output_NumSubnetworks++;
	}
	if (writeOutputFlag) {
		fprintf(stderr, "%llu subnetworks seeded by altered genes have been constructed and kept.\n", output_NumSubnetworks);
	}
	candidateSubnetworks.clear();
}

void constructSubnetworkCoresFromSet(const char * filename, SubnetworkEntry * & output_CandidateSubnetworks, llu & output_NumSubnetworks, Entry & patients, Entry & nodes, Entry & alterations, Bitmask *** geneSampleAlterations, bool writeOutputFlag = true) {
	int timerStart;
	vector<SubnetworkEntry> candidateSubnetworks;

	FILE * fin = fopen(filename, "r");
	char name[5000];
	char type[5000];
	char line[50000];
	char gene[1000];
	char alteration[1000];
	while (fscanf(fin, "%[^\t]\t%[^\t]\t", name, type) == 2) {
		if (writeOutputFlag) {
			fprintf(stderr, "%s\t", name);
		}
		bool doAnd = (strcmp(type, "I") == 0);	// Doing union by default
		fgets(line, 49999, fin);
		int readSoFar = 0;
		int chunkSize;
		map<int, llu> coreNodes;
		Bitmask coreSamples(patients.getSize());
		if (doAnd)
			coreSamples.negateSelf();
		while (sscanf(line + readSoFar, "%s%s%n", gene, alteration, &chunkSize) == 2) {
			readSoFar += chunkSize;
			if (writeOutputFlag) {
				fprintf(stderr, "[%s][%d][%s][%d]\t", gene, nodes.indices.count(string(gene)), alteration, alterations.indices.count(string(alteration)));
			}
			if (!nodes.indices.count(string(gene)) || !alterations.indices.count(string(alteration))) continue;
			int nodeIdx = nodes.indices[string(gene)];
			int altIdx = alterations.indices[string(alteration)];
			if (geneSampleAlterations[altIdx][nodeIdx] != nullptr) {
				llu mask = (llu(1) << altIdx);
				if (!coreNodes.count(nodeIdx))
					coreNodes[nodeIdx] = mask;
				else 
					coreNodes[nodeIdx] |= mask;
				if (doAnd)
					coreSamples.andWith(geneSampleAlterations[altIdx][nodeIdx]);
				else
					coreSamples.orWith(geneSampleAlterations[altIdx][nodeIdx]);
			}
		}
		if (writeOutputFlag) {
			fprintf(stderr, "\n");
		}
		if (!coreNodes.empty() && coreSamples.getSize() > 1) {
			SubnetworkEntry newEntry(coreNodes, coreSamples);
			newEntry.setName(name);
			candidateSubnetworks.push_back(newEntry);
		}
	}
	fclose(fin);
	if (writeOutputFlag) {
		fprintf(stderr, "Constructed core set subnetworks.\n");
	}
	
	output_CandidateSubnetworks = new SubnetworkEntry[candidateSubnetworks.size()];
	output_NumSubnetworks = 0;
	for (int i = 0; i < candidateSubnetworks.size(); i++) {
		output_CandidateSubnetworks[output_NumSubnetworks] = candidateSubnetworks[i];
		output_NumSubnetworks++;
	}
	if (writeOutputFlag) {
		fprintf(stderr, "%llu subnetworks seeded by altered genes have been constructed and kept.\n", output_NumSubnetworks);
	}
	candidateSubnetworks.clear();
}

bool cmpLessThan(int const & A, int const & B) {
	return A > B;
}

void updatePaths(const int startNode, const int startAlt, const map<int, llu> & nodesInSubnet, Bitmask *** overlap, Heap<int> & maxHeap, pair<int, int> * nodeStack, int * needsUpdating, Entry & patients, Entry & nodes, Entry & alterations, Bitmask *** geneSampleAlterations, Graph & G, unordered_map<string, Bitmask> & edgePatients) {
	// needsUpdating[idx] = 0 - doesn't need updating
	// needsUpdating[idx] = 1 - needs updating and it's already on stack
	// needsUpdating[idx] = 2 - needs updating and isn't on stack yet
	int stackSize = 0;
	nodeStack[stackSize++] = make_pair(startNode, startAlt);
	int startEntryIdx = startNode * alterations.getSize() + startAlt;
	needsUpdating[startEntryIdx] = 1;
	while (stackSize > 0 ) {
		int nodeIdx = nodeStack[stackSize - 1].first;
		int nodeAltIdx = nodeStack[stackSize - 1].second;
		stackSize--;
		int nodeEntryIdx = nodeIdx * alterations.getSize() + nodeAltIdx;
		if (needsUpdating[nodeEntryIdx] == 1) {	// Making sure that it wasn't already processed
			for (int j = 0; j < G.NSize[nodeIdx]; j++) {
				int neighbourIdx = G.N[nodeIdx][j];
				for (int neighbourAltIdx = 0; neighbourAltIdx < alterations.getSize(); neighbourAltIdx++) {
					if (geneSampleAlterations[neighbourAltIdx][neighbourIdx] != nullptr) {	// If it's even altered
						if (overlap[neighbourAltIdx][neighbourIdx] == nullptr)
							overlap[neighbourAltIdx][neighbourIdx] = new Bitmask(patients.getSize());
						string edgeName = nodes.names[nodeIdx] + " " + alterations.names[nodeAltIdx] + " " + nodes.names[neighbourIdx] + " " + alterations.names[neighbourAltIdx];
						if (edgePatients.count(edgeName)) {	// If this edge even exists with the chosen combination of node alterations
							Bitmask additionalPatients(edgePatients.find(edgeName)->second);	// Take the patients of the edge
							additionalPatients.andWith(overlap[nodeAltIdx][nodeIdx]);	// Intersect with the incoming node/alt's patients
							additionalPatients.andWith(geneSampleAlterations[neighbourAltIdx][neighbourIdx]);	// Intersect with patients that the destination node is altered in
							int oldNumPatientsWithPaths = overlap[neighbourAltIdx][neighbourIdx]->getSize();
							overlap[neighbourAltIdx][neighbourIdx]->orWith(additionalPatients);	// Add new paths

							// if (nodeIdx == nodes.indices["loss_3p"] && neighbourIdx == nodes.indices["VHL_3"] && nodeAltIdx == alterations.indices["CNLOSS"] && neighbourAltIdx == alterations.indices["INAC"]) {
							// 	fprintf(stderr, "\n\n");
							// 	fprintf(stderr, "%30s%d\n", "loss_3p overlap:", overlap[nodeAltIdx][nodeIdx]->getSize());
							// 	fprintf(stderr, "%30s%d\n", "Edge patients:", edgePatients.find(edgeName)->second.getSize());
							// 	fprintf(stderr, "%30s%d\n", "Additional patients:", additionalPatients.getSize());
							// 	fprintf(stderr, "%30s%d\n", "VHL overlap:", overlap[neighbourAltIdx][neighbourIdx]->getSize());
							// 	fprintf(stderr, "\n\n");
							// }
							// if (neighbourIdx == nodes.indices["VHL_3"] && neighbourAltIdx == alterations.indices["INAC"]) {
							// 	fprintf(stderr, "\n\n");
							// 	fprintf(stderr, "%30s %s %s\n", "Parent node:", nodes.names[nodeIdx].c_str(), alterations.names[nodeAltIdx].c_str());
							// 	fprintf(stderr, "%30s %d\n", "Parent overlap:", overlap[nodeAltIdx][nodeIdx]->getSize());
							// 	fprintf(stderr, "%30s %d\n", "Edge patients:", edgePatients.find(edgeName)->second.getSize());
							// 	fprintf(stderr, "%30s %d\n", "Additional patients:", additionalPatients.getSize());
							// 	fprintf(stderr, "%30s %d\n", "VHL overlap:", overlap[neighbourAltIdx][neighbourIdx]->getSize());
							// 	// fprintf(stderr, "%30s %s\n", "Contains R_K097?", overlap[neighbourAltIdx][neighbourIdx]->getBit(patients.indices["R_K097"]) == 1 ? "Yes." : "No.");
							// 	fprintf(stderr, "\n\n");
							// }
							if (oldNumPatientsWithPaths < overlap[neighbourAltIdx][neighbourIdx]->getSize()) {
								int entryIdx = neighbourIdx * alterations.getSize() + neighbourAltIdx;
								maxHeap.change_val(entryIdx, overlap[neighbourAltIdx][neighbourIdx]->getSize());
								if (nodesInSubnet.count(neighbourIdx) && needsUpdating[entryIdx] != 1)	// If the neighbour was inside the current subnetwork, then we need to update the rest of the subnetwork
									needsUpdating[entryIdx] = 2;	// But just mark it and we'll move into the node later when the rest of the neighbours have been updated
							}
						}
					}
				}
			}
			// After updating overlaps of all neighbours, now we recursively push the updates further
			for (int j = 0; j < G.NSize[nodeIdx]; j++) {
				int neighbourIdx = G.N[nodeIdx][j];
				for (int neighbourAltIdx = 0; neighbourAltIdx < alterations.getSize(); neighbourAltIdx++) {
					int entryIdx = neighbourIdx * alterations.getSize() + neighbourAltIdx;
					if (needsUpdating[entryIdx] == 2) {	// We want to have it on the stack only once
						nodeStack[stackSize++] = make_pair(neighbourIdx, neighbourAltIdx);
						needsUpdating[entryIdx] = 1;
					}
				}
			}
			needsUpdating[nodeEntryIdx] = 0;
		}
	}
}

void extendSubnetwork(const SubnetworkEntry & SN, int minDepth, SubnetworkEntry & outputSN, Entry & patients, Entry & nodes, Entry & alterations, Bitmask *** geneSampleAlterations, Graph & G, unordered_map<string, Bitmask> & edgePatients, Bitmask *** overlap = nullptr) {
	int timerStart;
	SubnetworkEntry newEntry(patients.getSize());
	const map<int, llu> & nodesInCore = SN.getNodesRef();
	const map<int, llu> & nodesInSubnet = newEntry.getNodesRef();
	bool shouldCleanOverlap = false;
	if (overlap == nullptr) {
		overlap = new Bitmask **[alterations.getSize()];
		for (int altIdx = 0; altIdx < alterations.getSize(); altIdx++) {
			overlap[altIdx] = new Bitmask * [nodes.getSize()];
			for (int nodeIdx = 0; nodeIdx < nodes.getSize(); nodeIdx++) {
				overlap[altIdx][nodeIdx] = nullptr;
			}
		}
		shouldCleanOverlap = true;
	}
	int * overlapSize = new int [nodes.getSize() * alterations.getSize()];
	pair<int, int> * nodeStack = new pair<int, int> [nodes.getSize() * alterations.getSize()];
	int * needsUpdating = new int [nodes.getSize() * alterations.getSize()];
	bool * isInSubnet = new bool [nodes.getSize() * alterations.getSize()];
	memset(overlapSize, 0, sizeof(overlapSize[0]) * (nodes.getSize() * alterations.getSize()));
	memset(needsUpdating, 0, sizeof(needsUpdating[0]) * (nodes.getSize() * alterations.getSize()));
	memset(isInSubnet, 0, sizeof(isInSubnet[0]) * (nodes.getSize() * alterations.getSize()));
	Heap<int> maxHeap(overlapSize, nodes.getSize() * alterations.getSize(), cmpLessThan);

	// Add core nodes one by one and update overlap sizes for all nodes in the subnetwork and outside neighbours
	for (auto it : nodesInCore) {
		int nodeIdx = it.first;
		Bitmask nodeAlterationMask(64);
		nodeAlterationMask.copylluBitmask(it.second);	// Go through colours of the node inside the subnetwork
		while (nodeAlterationMask.getSize()) {
			int altIdx = nodeAlterationMask.extractLowestOrderSetBitIndex();
			overlap[altIdx][nodeIdx] = new Bitmask(geneSampleAlterations[altIdx][nodeIdx]);
			// overlap[altIdx][nodeIdx]->andWith(SN.getSamplesPtr());
			newEntry.addNode(nodeIdx, altIdx);	// Add it to the subnetwork
			int heapIdx = nodeIdx * alterations.getSize() + altIdx;
			isInSubnet[heapIdx] = true;
			maxHeap.change_val(heapIdx, overlap[altIdx][nodeIdx]->getSize());	// And update the size of the overlap in the max heap
			updatePaths(nodeIdx, altIdx, nodesInSubnet, overlap, maxHeap, nodeStack, needsUpdating, patients, nodes, alterations, geneSampleAlterations, G, edgePatients);	// Update paths of all nodes in the subnet and direct neighbours outside
		}
	}
	// fprintf(stderr, "\n");
	// for (int g = 0; g < G.V; g++) {
	// 	for (int altIdx = 0; altIdx < alterations.getSize(); altIdx++) {
	// 		if (overlap[altIdx][g] != nullptr) {
	// 			// if (overlap[altIdx][g]->getSize() >= minDepth)
	// 				// fprintf(stderr, " %d", overlap[altIdx][g]->getSize());
	// 		}
	// 	}
	// }

	// Pick new neighbours as long as there are any with sufficient overlap
	while (overlapSize[maxHeap.getTopElementArrayIndex()] >= minDepth) {
		int heapIdx = maxHeap.getTopElementArrayIndex();
		int nodeIdx = heapIdx / alterations.getSize();
		int nodeAltIdx = heapIdx % alterations.getSize();
		maxHeap.pop();	// Remove the node/alt combo from the heap
		if (!isInSubnet[heapIdx]) {	// Core nodes are already in the subnet
			newEntry.addNode(nodeIdx, nodeAltIdx);	// Add it to the subnetwork
			updatePaths(nodeIdx, nodeAltIdx, nodesInSubnet, overlap, maxHeap, nodeStack, needsUpdating, patients, nodes, alterations, geneSampleAlterations, G, edgePatients);	// Update paths of all nodes in the subnet and direct neighbours outside
		}
	}

	outputSN = newEntry;
	// fprintf(stderr, "%d nodes in the subnetwork.\n", outputSN.getSize());

	// Deallocating memory
	if (shouldCleanOverlap) {
		for (int altIdx = 0; altIdx < alterations.getSize(); altIdx++) {
			for (int nodeIdx = 0; nodeIdx < G.V; nodeIdx++) {
				delete overlap[altIdx][nodeIdx];
			}
			delete overlap[altIdx];
		}
		delete overlap;
	}
	delete overlapSize;
	delete [] nodeStack;
	delete needsUpdating;
}

void constructSubnetworks(const SubnetworkEntry * coreSubnetworks, const llu numCores, const double minDepthRate, SubnetworkEntry * outputSubnetworks, Entry & patients, Entry & nodes, Entry & alterations, Bitmask *** geneSampleAlterations, Graph & G, unordered_map<string, Bitmask> & edgePatients, bool writeOutputFlag = true) {
	int timerStart;
	// fprintf(stderr, "Extending core subnetworks.\n");
	timerStart = clock();
	for (int i = 0, lastProg = 0; i < numCores; i++) {
		if (writeOutputFlag) {
			int progress = 1000 * double(i + 1) / double(numCores);
			if (progress > lastProg) {
				fprintf(stderr, "\r\t%.1lf%%", double(progress)/10.0);
				lastProg = progress;
			}
		}
		// extendSubnetwork(coreSubnetworks[i], minDepthRate * coreSubnetworks[i].getNumSamples(), outputSubnetworks[i], patients, nodes, alterations, geneSampleAlterations, G, edgePatients);
		extendSubnetwork(coreSubnetworks[i], minDepthRate, outputSubnetworks[i], patients, nodes, alterations, geneSampleAlterations, G, edgePatients);
	}
	if (writeOutputFlag) {
		fprintf(stderr, "\r\tDone. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC);
	}
}

int chooseOptimalSubnetwork(const SubnetworkEntry * subnetworks, const llu numSubnets, const double minDepthRate, const double edgeThreshold, const SubnetworkEntry * cores, string outFolder, Entry & patients, Entry & nodes, Entry & alterations, Bitmask *** geneSampleAlterations, Graph & G, unordered_map<string, Bitmask> & edgePatients) {
	// string outDistributionRecurrence	= outFolder + "/subnetwork_recurrence.txt";
	string outSubnFolder				= outFolder + "/subnetworks";
	char command[1000];
	sprintf(command, "mkdir -p %s", outSubnFolder.c_str());
	system(command);
	// Prepare the overlap data structure
	Bitmask *** overlap = new Bitmask **[alterations.getSize()];
	for (int altIdx = 0; altIdx < alterations.getSize(); altIdx++) {
		overlap[altIdx] = new Bitmask * [G.V];
		for (int nodeIdx = 0; nodeIdx < G.V; nodeIdx++) {
			overlap[altIdx][nodeIdx] = new Bitmask(patients.getSize());
		}
	}

	FILE * foutSizes = fopen((outFolder + "/subnetSizes.txt").c_str(), "w");
	fprintf(foutSizes, "Core\tSize\n");
	int maxIdx = 0;
	int maxCnt = 1;
	for (int i = 1; i < numSubnets; i++) {
		// Print the core nodes
		const map<int, llu> & nodesInCore = cores[i].getNodesRef();
		if (cores[i].getName().size() > 0) {
			fprintf(foutSizes, "%s\t", cores[i].getName().c_str());
		}
		else {
			for (auto it : nodesInCore) {
				int nodeIdx = it.first;
				llu altMask = it.second;
				fprintf(foutSizes, "%s", nodes.names[nodeIdx].c_str());
				Bitmask alterationMask(64);
				alterationMask.copylluBitmask(altMask);	// Go through alterations of the node
				while (alterationMask.getSize()) {
					int altIdx = alterationMask.extractLowestOrderSetBitIndex();
					fprintf(foutSizes, "|%s", alterations.names[altIdx].c_str());
				}
				fprintf(foutSizes, " ");
			}
		}
		fprintf(foutSizes, "%d\n", subnetworks[i].getSize());
		if (subnetworks[i].getSize() > subnetworks[maxIdx].getSize()) {
			maxIdx = i;
			maxCnt = 1;
		}
		else if (subnetworks[i].getSize() == subnetworks[maxIdx].getSize()) {
			maxCnt++;
		}
	}
	int maxSize = subnetworks[maxIdx].getSize();
	fclose(foutSizes);
	fprintf(stderr, "Maximum subnetwork size is %d.\n", maxSize);
	fprintf(stderr, "There are %d subnetworks of maximum size.\n", maxCnt);
	fprintf(stderr, "The first such subnetwork has index %d.\n", maxIdx);

	// Printing maximum-sized subnetworks
	int printIdx = 0;
	for (int subnetIdx = 0; subnetIdx < numSubnets; subnetIdx++) {
		if (subnetworks[subnetIdx].getSize() == maxSize) {
			printIdx++;
			char filename[1000];
			sprintf(filename, "%s/%d.txt", outSubnFolder.c_str(), printIdx);
			FILE * foutSubnet = fopen(filename, "w");
			sprintf(filename, "%s/%d_nodes.txt", outSubnFolder.c_str(), printIdx);
			FILE * foutNodes = fopen(filename, "w");
			sprintf(filename, "%s/%d_edges.txt", outSubnFolder.c_str(), printIdx);
			FILE * foutEdges = fopen(filename, "w");
			sprintf(filename, "%s/%d_adj.txt", outSubnFolder.c_str(), printIdx);
			FILE * foutAdj = fopen(filename, "w");
			sprintf(filename, "%s/%d_core.txt", outSubnFolder.c_str(), printIdx);
			FILE * foutCore = fopen(filename, "w");
			sprintf(filename, "%s/%d_patients.txt", outSubnFolder.c_str(), printIdx);
			FILE * foutCoreSamples = fopen(filename, "w");
			sprintf(filename, "%s/%d_patientMap.txt", outSubnFolder.c_str(), printIdx);
			FILE * foutPatientMap = fopen(filename, "w");
			sprintf(filename, "%s/%d_tree_adj.txt", outSubnFolder.c_str(), printIdx);
			FILE * foutTree = fopen(filename, "w");

			// Print the core nodes
			const map<int, llu> & nodesInCore = cores[subnetIdx].getNodesRef();
			if (cores[subnetIdx].getName().size() > 0)
				fprintf(foutSubnet, "%s\n", cores[subnetIdx].getName().c_str());
			fprintf(foutSubnet, "%d\n", (int)nodesInCore.size());
			for (auto it : nodesInCore) {
				int nodeIdx = it.first;
				llu altMask = it.second;
				fprintf(foutSubnet, "%s", nodes.names[nodeIdx].c_str());
				fprintf(foutCore, "%s", nodes.names[nodeIdx].c_str());
				Bitmask alterationMask(64);
				alterationMask.copylluBitmask(altMask);	// Go through alterations of the node
				while (alterationMask.getSize()) {
					int altIdx = alterationMask.extractLowestOrderSetBitIndex();
					fprintf(foutSubnet, "\t%s", alterations.names[altIdx].c_str());
					fprintf(foutCore, "\t%s", alterations.names[altIdx].c_str());
				}
				fprintf(foutSubnet, "\n");
				fprintf(foutCore, "\n");
			}

			// Print the core patients
			Bitmask patientsInCore(cores[subnetIdx].getSamplesPtr());
			fprintf(foutSubnet, "%d\n", patientsInCore.getSize());
			while (patientsInCore.getSize()) {
				int sampleIdx = patientsInCore.extractLowestOrderSetBitIndex();
				fprintf(foutCoreSamples, "%s\n", patients.names[sampleIdx].c_str());
				fprintf(foutSubnet, "%s ", patients.names[sampleIdx].c_str());
			}
			fprintf(foutSubnet, "\n");

			// Find and print nodes and their alterations within the subnetwork
			const map<int, llu> & nodesInSubnet = subnetworks[subnetIdx].getNodesRef();
			vector<int> nodeIdxInSubnet(nodesInSubnet.size(), -1);
			unordered_map<int, int> adjIdx;
			fprintf(foutSubnet, "%d\n", (int)nodesInSubnet.size());
			fprintf(stderr, "> | %d | ", cores[subnetIdx].getSamplesPtr()->getSize());
			for (auto it : nodesInSubnet) {
				int nodeIdx = it.first;
				llu altMask = it.second;
				int currentIdx = adjIdx.size();
				adjIdx[nodeIdx] = currentIdx;
				nodeIdxInSubnet[currentIdx] = nodeIdx;

				fprintf(foutNodes, "%s", nodes.names[nodeIdx].c_str());
				fprintf(stderr, "%s", nodes.names[nodeIdx].c_str());
				fprintf(foutSubnet, "%s", nodes.names[nodeIdx].c_str());
				// if (geneChrLoc.count(nodes.names[nodeIdx]))
					// fprintf(foutSubnet, "\t%s", geneChrLoc[nodes.names[nodeIdx]].c_str());
				Bitmask alterationMask(64);
				alterationMask.copylluBitmask(altMask);	// Go through alterations of the node
				while (alterationMask.getSize()) {
					int altIdx = alterationMask.extractLowestOrderSetBitIndex();
					fprintf(foutNodes, "\t%s", alterations.names[altIdx].c_str());
					fprintf(stderr, " (%s)\t", alterations.names[altIdx].c_str());
					fprintf(foutSubnet, "\t%s", alterations.names[altIdx].c_str());
					// overlap[altIdx][nodeIdx]->clear();
				}
				fprintf(foutNodes, "\n");
				fprintf(foutSubnet, "\n");
			}
			fprintf(stderr, "\n");

			// Recalculate covered patients of the subnetwork by recomputing it
			SubnetworkEntry tempNet;
			for (int altIdx = 0; altIdx < alterations.getSize(); altIdx++) {
				for (int nodeIdx = 0; nodeIdx < nodes.getSize(); nodeIdx++) {
					overlap[altIdx][nodeIdx]->clear();
				}
			}
			// extendSubnetwork(cores[subnetIdx], minDepthRate * cores[subnetIdx].getNumSamples(), tempNet, patients, nodes, alterations, geneSampleAlterations, G, edgePatients, overlap);
			extendSubnetwork(cores[subnetIdx], minDepthRate, tempNet, patients, nodes, alterations, geneSampleAlterations, G, edgePatients, overlap);
			// for (auto it : nodesInSubnet) {
			// 	int nodeIdx = it.first;
			// 	llu altMask = it.second;
			// 	Bitmask alterationMask(64);
			// 	alterationMask.copylluBitmask(altMask);	// Go through alterations of the node
			// 	while (alterationMask.getSize()) {
			// 		int altIdx = alterationMask.extractLowestOrderSetBitIndex();
			// 		fprintf(stderr, "%s\t%s\t%d\n", nodes.names[nodeIdx].c_str(), alterations.names[altIdx].c_str(), overlap[altIdx][nodeIdx]->getSize());
			// 	}
			// }
			// fprintf(stderr, "%.2lf\t%d\t%d\n", minDepthRate, cores[subnetIdx].getNumSamples(), int(minDepthRate * cores[subnetIdx].getNumSamples()));
			// for (int altIdx = 0; altIdx < alterations.getSize(); altIdx++) {
			// 	for (int nodeIdx = 0; nodeIdx < nodes.getSize(); nodeIdx++) {
			// 		if (overlap[altIdx][nodeIdx]->getSize() < int(minDepthRate * cores[subnetIdx].getNumSamples()))
			// 			// fprintf(stderr, "%s %s\n", nodes.names[nodeIdx].c_str(), alterations.names[altIdx].c_str());
			// 			overlap[altIdx][nodeIdx]->clear();
			// 	}
			// }
			// fprintf(stderr, "\n\n");
			// int VHLIdx = nodes.indices["VHL_3"];
			// int inacIdx = alterations.indices["INAC"];
			// int lap = overlap[inacIdx][VHLIdx] -> getSize();
			// fprintf(stderr, "VHL overlap to loss_3p is %d\n", lap);
			// fprintf(stderr, "\n\n");
			// for (auto it : nodesInSubnet) {
			// 	int nodeIdx = it.first;
			// 	llu altMask = it.second;
			// 	Bitmask alterationMask(64);
			// 	alterationMask.copylluBitmask(altMask);	// Go through alterations of the node
			// 	while (alterationMask.getSize()) {
			// 		int altIdx = alterationMask.extractLowestOrderSetBitIndex();
			// 		fprintf(stderr, "%s\t%s\t%d\n", nodes.names[nodeIdx].c_str(), alterations.names[altIdx].c_str(), overlap[altIdx][nodeIdx]->getSize());
			// 	}
			// }
			// fprintf(stderr, "\n");

			// Detect edges
			double ** adj = new double * [nodesInSubnet.size()];
			double ** W = new double * [nodesInSubnet.size()];
			double ** Agreed = new double * [nodesInSubnet.size()];
			Bitmask *** SubnetworkPatientsWithEdge = new Bitmask ** [nodesInSubnet.size()];
			Bitmask ** ownPatients = new Bitmask * [nodesInSubnet.size()];
			for (int nodeIdx = 0; nodeIdx < nodesInSubnet.size(); nodeIdx++) {
				SubnetworkPatientsWithEdge[nodeIdx] = new Bitmask *[nodesInSubnet.size()];
				for (int neighbourIdx = 0; neighbourIdx < nodesInSubnet.size(); neighbourIdx++) {
					SubnetworkPatientsWithEdge[nodeIdx][neighbourIdx] = new Bitmask(patients.getSize());
				}
				ownPatients[nodeIdx] = new Bitmask(patients.getSize());
			}
			for (auto it : adjIdx) {
				int idxInMatrix = it.second;
				adj[idxInMatrix] = new double [nodesInSubnet.size()];
				W[idxInMatrix] = new double [nodesInSubnet.size()];
				Agreed[idxInMatrix] = new double [nodesInSubnet.size()];
				memset(adj[idxInMatrix], 0, sizeof(adj[idxInMatrix][0]) * nodesInSubnet.size());
				memset(W[idxInMatrix], 0, sizeof(W[idxInMatrix][0]) * nodesInSubnet.size());
				memset(Agreed[idxInMatrix], 0, sizeof(Agreed[idxInMatrix][0]) * nodesInSubnet.size());
			}

			// Patient graphs
			patientsInCore = cores[subnetIdx].getSamplesPtr();
			while (patientsInCore.getSize()) {
				int sampleIdx = patientsInCore.extractLowestOrderSetBitIndex();
			// for (int sampleIdx = 0; sampleIdx < patients.getSize(); sampleIdx++) {
				bool thereIsAnEdge = false;
				for (int i = 0; i < nodesInSubnet.size(); i++) {
					memset(adj[i], 0, sizeof(adj[i][0]) * nodesInSubnet.size());
				}
				for (auto it : nodesInSubnet) {
					int nodeIdx = it.first;
					int idxInMatrix1 = adjIdx[nodeIdx];
					llu altMask = it.second;
					Bitmask alterationMask(64);
					alterationMask.copylluBitmask(altMask);	// Go through alterations of the node
					while (alterationMask.getSize()) {
						int altIdx = alterationMask.extractLowestOrderSetBitIndex();
						if (overlap[altIdx][nodeIdx]->getBit(sampleIdx)) {
							for (int j = 0; j < G.NSize[nodeIdx]; j++) {
								int neighbourIdx = G.N[nodeIdx][j];
								if (nodesInSubnet.count(neighbourIdx)) {
									int idxInMatrix2 = adjIdx[neighbourIdx];
									for (int neighbourAltIdx = 0; neighbourAltIdx < alterations.getSize(); neighbourAltIdx++) {
										string edgeName = nodes.names[nodeIdx] + " " + alterations.names[altIdx] + " " + nodes.names[neighbourIdx] + " " + alterations.names[neighbourAltIdx];
										if (edgePatients.count(edgeName)) {	// If this edge even exists with the chosen combination of node alterations
											Bitmask validPatients(edgePatients.find(edgeName)->second);	// Take the patients of the edge
											if (overlap[neighbourAltIdx][neighbourIdx]->getBit(sampleIdx) && validPatients.getBit(sampleIdx)) {
												adj[idxInMatrix1][idxInMatrix2] = 1;
												thereIsAnEdge = true;
											}
										}
									}
								}
							}
						}
					}
				}
				if (thereIsAnEdge) {
					fprintf(foutAdj, "\n%s\n", patients.names[sampleIdx].c_str());
					for (int i = 0; i < nodesInSubnet.size(); i++) {
						for (int j = 0; j < nodesInSubnet.size(); j++) {
							fprintf(foutAdj, "%d ", (int)adj[i][j]);
						}
						fprintf(foutAdj, "\n");
					}
				}
			}

			// Union graph
			for (int i = 0; i < nodesInSubnet.size(); i++) {
				memset(adj[i], 0, sizeof(adj[i][0]) * nodesInSubnet.size());
			}
			for (auto it : nodesInSubnet) {
				int nodeIdx = it.first;
				for (int j = 0; j < G.NSize[nodeIdx]; j++) {
					int neighbourIdx = G.N[nodeIdx][j];
					if (nodesInSubnet.count(neighbourIdx)) {	// We only care about subnetwork nodes
						// adj[ adjIdx[nodeIdx] ][ adjIdx[neighbourIdx] ] = 1;
						fprintf(foutEdges, "%s\t%s\n", nodes.names[nodeIdx].c_str(), nodes.names[neighbourIdx].c_str());	// Print the edge
						// int totalCnt = 0;
						// int presentIn = 0;
						llu altMask = it.second;
						Bitmask alterationMask(64);
						alterationMask.copylluBitmask(altMask);	// Go through alterations of the node
						while (alterationMask.getSize()) {
							int altIdx = alterationMask.extractLowestOrderSetBitIndex();
							llu neighbourAltMask = nodesInSubnet.find(neighbourIdx)->second;
							Bitmask neighbourAlterationMask(64);
							neighbourAlterationMask.copylluBitmask(neighbourAltMask);	// Go through alterations of the neighbour
							while (neighbourAlterationMask.getSize()) {
								int neighbourAltIdx = neighbourAlterationMask.extractLowestOrderSetBitIndex();
								string edgeName = nodes.names[nodeIdx] + " " + alterations.names[altIdx] + " " + nodes.names[neighbourIdx] + " " + alterations.names[neighbourAltIdx];
								if (edgePatients.count(edgeName)) {	// If this edge even exists with the chosen combination of node alterations
									Bitmask validPatients(edgePatients.find(edgeName)->second);	// Take the patients of the edge
									// patientsInCore = cores[subnetIdx].getSamplesPtr();
									// while (patientsInCore.getSize()) {
									// 	int sampleIdx = patientsInCore.extractLowestOrderSetBitIndex();
									// 	if (validPatients.getBit(sampleIdx) && overlap[altIdx][nodeIdx]->getBit(sampleIdx) && overlap[neighbourAltIdx][neighbourIdx]->getBit(sampleIdx)) {
									// 		adj[ adjIdx[nodeIdx] ][ adjIdx[neighbourIdx] ]++;
									// 		// SubnetworkPatientsWithEdge[ adjIdx[nodeIdx] ][ adjIdx[neighbourIdx] ] -> setBit(sampleIdx, 1);
									// 		// SubnetworkPatientsWithEdge[ adjIdx[neighbourIdx] ][ adjIdx[nodeIdx] ] -> setBit(sampleIdx, 1);
									// 		Bitmask tempMask(overlap[altIdx][nodeIdx]);
									// 		tempMask.andWith(overlap[neighbourAltIdx][neighbourIdx]);
									// 		SubnetworkPatientsWithEdge[ adjIdx[nodeIdx] ][ adjIdx[neighbourIdx] ] -> orWith(tempMask);
									// 		// SubnetworkPatientsWithEdge[ adjIdx[neighbourIdx] ][ adjIdx[nodeIdx] ] -> orWith(tempMask);
									// 		ownPatients[ adjIdx[nodeIdx] ] -> orWith(overlap[altIdx][nodeIdx]);
									// 	}
									// }
									// validPatients.andWith(patientsInCore);	// redundant since overlap only contains patients in core
									validPatients.andWith(overlap[altIdx][nodeIdx]);
									validPatients.andWith(overlap[neighbourAltIdx][neighbourIdx]);
									adj[ adjIdx[nodeIdx] ][ adjIdx[neighbourIdx] ] += validPatients.getSize();
									SubnetworkPatientsWithEdge[ adjIdx[nodeIdx] ][ adjIdx[neighbourIdx] ] -> orWith(validPatients);
								}
							}
						}
					}
				}
				for (int altIdx = 0; altIdx < alterations.getSize(); altIdx++) {
					ownPatients[ adjIdx[nodeIdx] ] -> orWith(overlap[altIdx][nodeIdx]);
					// ownPatients[ adjIdx[nodeIdx] ] -> orWith(geneSampleAlterations[altIdx][nodeIdx]);
				}
			}

			for (int i = 0; i < nodesInSubnet.size(); i++) {
				for (int j = 0; j < nodesInSubnet.size(); j++) {
					if (SubnetworkPatientsWithEdge[i][j]->getSize()) {
						// Bitmask intersectionPatients(ownPatients[i]);
						// intersectionPatients.andWith(ownPatients[j]);
						// // Bitmask unionPatients(ownPatients[i]);
						// // unionPatients.orWith(ownPatients[j]);
						// Bitmask sameClone(SubnetworkPatientsWithEdge[i][j]);
						// sameClone.andWith(SubnetworkPatientsWithEdge[j][i]);
						// Bitmask reverseEdgesNotSameClone(SubnetworkPatientsWithEdge[j][i]);
						// reverseEdgesNotSameClone.xorWith(sameClone);

						// double Pij = double(SubnetworkPatientsWithEdge[i][j]->getSize()) / double(intersectionPatients.getSize());
						// double cij = double(reverseEdgesNotSameClone.getSize()) / double(intersectionPatients.getSize());
						// double Aij = adj[i][j] / double(ownPatients[j] -> getSize());
						// // W[i][j] = Pij - cij - Aij;
						// W[i][j] = Pij - cij;
						// Agreed[i][j] = Pij;
						// if (Aij < 0.3) {
						// 	W[i][j] = 0;
						// }

						// Bitmask unionPatients(ownPatients[i]);
						// unionPatients.orWith(ownPatients[j]);
						Bitmask intersectionPatients(ownPatients[i]);
						intersectionPatients.andWith(ownPatients[j]);
						// Bitmask separatePatients(unionPatients);
						// separatePatients.xorWith(intersectionPatients);

						// Bitmask firstAlone(ownPatients[i]);
						// firstAlone.xorWith(intersectionPatients);
						// Bitmask secondAlone(ownPatients[j]);
						// secondAlone.xorWith(intersectionPatients);

						Bitmask sameClone(SubnetworkPatientsWithEdge[i][j]);
						sameClone.andWith(SubnetworkPatientsWithEdge[j][i]);
						Bitmask reverseEdgesNotSameClone(SubnetworkPatientsWithEdge[j][i]);
						reverseEdgesNotSameClone.xorWith(sameClone);

						double Pij = double(SubnetworkPatientsWithEdge[i][j]->getSize()) / double(intersectionPatients.getSize());
						double cij = double(reverseEdgesNotSameClone.getSize()) / double(intersectionPatients.getSize());
						// double sij = double(separatePatients.getSize()) / double(unionPatients.getSize());
						// double Iij = double(secondAlone.getSize()) / double(ownPatients[j]->getSize());
						// double Aij = adj[i][j] / double(ownPatients[j] -> getSize());
						// W[i][j] = Pij - cij - Aij;
						// W[i][j] = Pij - cij - sij;
						// W[i][j] = Pij - cij - Iij;
						W[i][j] = Pij - cij;

						Agreed[i][j] = Pij;
						// if (Aij < 0.5) {
						// 	W[i][j] = 0;
						// }
					}
				}
			}
			// {
			// 	fprintf(stderr, "\n\n");
			// 	int loss3pIdxInMatrix = adjIdx[ nodes.indices["loss_3p"] ];
			// 	int VHLIdx = nodes.indices["VHL_3"];
			// 	int VHLIdxInMatrix = adjIdx[ nodes.indices["VHL_3"] ];
			// 	int inacIdx = alterations.indices["INAC"];
			// 	string edgeName = "loss_3p CNLOSS VHL_3 INAC";
			// 	int edgePatientCount = edgePatients.find(edgeName)->second.getSize();

			// 	Bitmask tempMask2(SubnetworkPatientsWithEdge[loss3pIdxInMatrix][VHLIdxInMatrix]);
			// 	Bitmask tempMask3(overlap[inacIdx][VHLIdx]);
			// 	fprintf(stderr, "\ncommon: %d\toverlap: %d\tadj: %d\tEdgePatientCount: %d\n", tempMask2.getSize(), tempMask3.getSize(), (int)adj[loss3pIdxInMatrix][VHLIdxInMatrix], edgePatientCount);
			// 	fprintf(stderr, "\n\n");
			// }
			// tempMask2.xorWith(tempMask3);
			// while (tempMask2.getSize()) {
			// 	int sampleIdx = tempMask2.extractLowestOrderSetBitIndex();
			// 	fprintf(stderr, "NO EDGE IN PATIENT %s\n", patients.names[sampleIdx].c_str());
			// }

			// Print the adjacency matrix
			// Union graph
			fprintf(foutAdj, "\n");
			for (int i = 0; i < nodesInSubnet.size(); i++) {
				for (int j = 0; j < nodesInSubnet.size(); j++) {
					// fprintf(foutAdj, "%.2lf ", adj[i][j]);
					fprintf(foutAdj, "%.0lf ", adj[i][j]);
				}
				fprintf(foutAdj, "\n");
			}
			fprintf(foutAdj, "\n");
			for (int i = 0; i < nodesInSubnet.size(); i++) {
				for (int j = 0; j < nodesInSubnet.size(); j++) {
					fprintf(foutAdj, "%.2lf ", W[i][j]);
				}
				fprintf(foutAdj, "\n");
			}
			fprintf(foutAdj, "\n");

			patientsInCore = cores[subnetIdx].getSamplesPtr();
			fprintf(foutPatientMap, "G/S");
			while (patientsInCore.getSize()) {
				int sampleIdx = patientsInCore.extractLowestOrderSetBitIndex();
				fprintf(foutPatientMap, "\t%s", patients.names[sampleIdx].c_str());
			}
			// Print the covered patients for every node
			fprintf(foutSubnet, "\n");
			for (auto it : nodesInSubnet) {
				unordered_set<int> nodePatients;
				int nodeIdx = it.first;
				// if (nodesInCore.count(nodeIdx))
				// 	continue;
				llu altMask = it.second;
				Bitmask alterationMask(64);
				alterationMask.copylluBitmask(altMask);	// Go through alterations of the node
				while (alterationMask.getSize()) {
					int altIdx = alterationMask.extractLowestOrderSetBitIndex();
					fprintf(foutSubnet, "%s\t%s\t%d", nodes.names[nodeIdx].c_str(), alterations.names[altIdx].c_str(), overlap[altIdx][nodeIdx]->getSize());
					Bitmask tempMask(overlap[altIdx][nodeIdx]);
					while (tempMask.getSize()) {	// Go through patients covered by this node with this alteration
						int sampleIdx = tempMask.extractLowestOrderSetBitIndex();
						fprintf(foutSubnet, "\t%s", patients.names[sampleIdx].c_str());
						nodePatients.insert(sampleIdx);
					}
					fprintf(foutSubnet, "\n");
				}

				fprintf(foutPatientMap, "\n%s", nodes.names[nodeIdx].c_str());
				patientsInCore = cores[subnetIdx].getSamplesPtr();
				while (patientsInCore.getSize()) {
					int sampleIdx = patientsInCore.extractLowestOrderSetBitIndex();
					if (nodePatients.count(sampleIdx)) {
						fprintf(foutPatientMap, "\t1");
					}
					else {
						fprintf(foutPatientMap, "\t0");
					}
				}
			}

			fclose(foutSubnet);
			fclose(foutNodes);
			fclose(foutEdges);
			fclose(foutAdj);
			fclose(foutCore);
			fclose(foutCoreSamples);
			fclose(foutPatientMap);
			
			//
			//	TREE CONSTRUCTION via ILP
			//
			
			GRBEnv * env			= 0;
			GRBLinExpr objective	= 0;
			GRBVar ** X				= 0;
			GRBVar ** PTH			= 0;
			GRBVar *** PTHA			= 0;
			GRBVar ** F				= 0;
			GRBVar ** PD			= 0;
			GRBVar ** A				= 0;
			GRBVar *** PA			= 0;
			GRBVar * D				= 0;
			// GRBVar * I				= 0;
			// double * ub = new double [nodesInSubnet.size()];

			int sourceIdx = 0;
			for (auto it : nodesInCore) {
				sourceIdx = adjIdx[it.first];
				break;
			}

			try {
				env = new GRBEnv("ILPlog.txt");
				GRBModel model = GRBModel(*env);

				X = new GRBVar * [nodesInSubnet.size()];
				PTH = new GRBVar * [nodesInSubnet.size()];
				PTHA = new GRBVar ** [nodesInSubnet.size()];
				A = new GRBVar * [nodesInSubnet.size()];
				PD = new GRBVar * [nodesInSubnet.size()];
				PA = new GRBVar ** [nodesInSubnet.size()];
				F = new GRBVar * [nodesInSubnet.size()];
				for (int i = 0; i < nodesInSubnet.size(); i++) {
					X[i] = model.addVars(nodesInSubnet.size(), GRB_BINARY);
					PTH[i] = model.addVars(patients.getSize(), GRB_CONTINUOUS);
					PTHA[i] = new GRBVar * [nodesInSubnet.size()];
					A[i] = model.addVars(nodesInSubnet.size(), GRB_CONTINUOUS);
					PD[i] = model.addVars(nodesInSubnet.size(), GRB_CONTINUOUS);
					PA[i] = new GRBVar * [nodesInSubnet.size()];
					F[i] = model.addVars(nodesInSubnet.size(), GRB_CONTINUOUS);
					for (int j = 0; j < nodesInSubnet.size(); j++) {
						PA[i][j] = model.addVars(nodesInSubnet.size(), GRB_CONTINUOUS);
						PTHA[i][j] = model.addVars(patients.getSize(), GRB_CONTINUOUS);
					}
				}
				D = model.addVars(nodesInSubnet.size(), GRB_CONTINUOUS);
				// I = model.addVars(nodesInSubnet.size(), GRB_CONTINUOUS);
				// for (int i = 0; i < nodesInSubnet.size(); i++) {	// Naming variables //////////////
				// 	char varName[100];
				// 	sprintf(varName, "D%d", i);
				// 	D[i] = model.addVar(0, nodesInSubnet.size() - 1, 0, GRB_CONTINUOUS, varName);
				// }
				// for (int i = 0; i < nodesInSubnet.size(); i++) {	// Naming variables //////////////
				// 	for (int j = 0; j < nodesInSubnet.size(); j++) {
				// 		if (W[i][j] >= edgeThreshold) {
				// 			char varName[100];
				// 			sprintf(varName, "X%d_%d", i, j);
				// 			X[i][j] = model.addVar(0, 1, 0, GRB_BINARY, varName);
				// 			sprintf(varName, "A%d_%d", i, j);
				// 			A[i][j] = model.addVar(0, 1, 0, GRB_BINARY, varName);
				// 			sprintf(varName, "F%d_%d", i, j);
				// 			F[i][j] = model.addVar(0, nodesInSubnet.size() - 1, 0, GRB_CONTINUOUS, varName);
				// 			sprintf(varName, "PD%d_%d", i, j);
				// 			PD[i][j] = model.addVar(0, nodesInSubnet.size() - 1, 0, GRB_CONTINUOUS, varName);
				// 		}
				// 	}
				// }

				// Objective function //
				for (int i = 0; i < nodesInSubnet.size(); i++) {
					// double w = 1.0 / double(nodesInSubnet.size() * nodesInSubnet.size());
					// objective += w * D[i];
					objective += D[i];
					// objective += I[i];
				}
				// for (int i = 0; i < nodesInSubnet.size(); i++) {
				// 	for (int j = 0; j < nodesInSubnet.size(); j++) {
				// 		if (adj[i][j] > 0) {
				// 			// double w = 3.0;
				// 			// objective += w * C[i][j] * X[i][j];
				// 			// double w = (nodesInSubnet.size() - 1) * (W[i][j] - W[j][i]);
				// 			// double w = (nodesInSubnet.size()) * (C[i][j] - C[j][i]);
				// 			// objective += w * X[i][j];
				// 			objective -= A[i][j];
				// 		}
				// 	}
				// }
				model.setObjective(objective, GRB_MAXIMIZE);
				fprintf(stderr, "\tConstructed the objective function.\n");

				// Constraint 1 //
				model.addConstr(D[sourceIdx] == 0);

				// Constraint 1.5 //
				for (int j = 0; j < nodesInSubnet.size(); j++) {
					model.addConstr(X[j][sourceIdx] == 0);
					// model.addConstr(F[j][sourceIdx] == 0);
				}

				// Constraint 2 //
				for (int i = 0; i < nodesInSubnet.size(); i++) {
					for (int j = 0; j < nodesInSubnet.size(); j++) {
						if (j != sourceIdx) {
							model.addConstr(PD[i][j] <= D[i]);
						}
					}
				}

				// Constraint 3 //
				for (int i = 0; i < nodesInSubnet.size(); i++) {
					for (int j = 0; j < nodesInSubnet.size(); j++) {
						if (j != sourceIdx) {
							model.addConstr(PD[i][j] <= nodesInSubnet.size() * X[i][j]);
						}
					}
				}

				// Constraint 4 //
				for (int j = 0; j < nodesInSubnet.size(); j++) {
					if (j != sourceIdx) {
						GRBLinExpr expr = D[j];
						for (int i = 0; i < nodesInSubnet.size(); i++) {
							expr -= PD[i][j];
						}
						model.addConstr(expr == 1);
					}
				}

				// Constraint 5 //
				for (int j = 0; j < nodesInSubnet.size(); j++) {
					if (j != sourceIdx) {
						GRBLinExpr expr = 0;
						for (int i = 0; i < nodesInSubnet.size(); i++) {
							expr += X[i][j];
						}
						model.addConstr(expr == 1);
					}
				}

				// Constraint 6 //
				{
					GRBLinExpr expr = 0;
					for (int j = 0; j < nodesInSubnet.size(); j++) {
						expr += F[sourceIdx][j];
					}
					model.addConstr(expr == nodesInSubnet.size() - 1);
				}

				// Constraint 7 //
				for (int j = 0; j < nodesInSubnet.size(); j++) {
					if (j != sourceIdx) {
						GRBLinExpr expr = 0;
						for (int i = 0; i < nodesInSubnet.size(); i++) {
							expr += F[i][j];
						}
						for (int i = 0; i < nodesInSubnet.size(); i++) {
							expr -= F[j][i];
						}
						model.addConstr(expr == 1);
					}
				}

				// Constraint 8 //
				for (int i = 0; i < nodesInSubnet.size(); i++) {
					for (int j = 0; j < nodesInSubnet.size(); j++) {
						model.addConstr(nodesInSubnet.size() * X[i][j] >= F[i][j]);
					}
				}

				// Constraint 9 //
				for (int i = 0; i < nodesInSubnet.size(); i++) {
					for (int j = 0; j < nodesInSubnet.size(); j++) {
						model.addConstr(X[i][j] <= F[i][j]);
					}
				}

				// Constraint 10 //
				for (int i = 0; i < nodesInSubnet.size(); i++) {
					for (int j = 0; j < nodesInSubnet.size(); j++) {
						if (W[i][j] < edgeThreshold) {
						// if (W[i][j] < edgeThreshold || adj[i][j] < 0.3 * adj[sourceIdx][j]) {
							// double w = double(1) / double(nodesInSubnet.size());
							model.addConstr(X[i][j] == 0);
							model.addConstr(A[i][j] == 0);
						}
					}
				}

				// Constraint 11 //
				for (int i = 0; i < nodesInSubnet.size(); i++) {
					for (int j = 0; j < nodesInSubnet.size(); j++) {
						model.addConstr(A[i][j] >= X[i][j]);
						model.addConstr(A[i][j] <= 1);

						GRBLinExpr expr = 0;
						for (int k = 0; k < nodesInSubnet.size(); k++) {
							// if (W[i][k] >= edgeThreshold && W[k][j] >= edgeThreshold) {
								model.addConstr(PA[i][j][k] <= A[i][k]);
								model.addConstr(PA[i][j][k] <= A[k][j]);
								model.addConstr(PA[i][j][k] >= A[i][k] + A[k][j] - 1);
								expr += PA[i][j][k];
							// }
						}
						model.addConstr(nodesInSubnet.size() * A[i][j] >= expr);
						model.addConstr(A[i][j] <= expr + X[i][j]);
					}
				}

				// // Constraint 12 //
				// for (int i = 0; i < nodesInSubnet.size(); i++) {
				// 	GRBLinExpr expr = 0;
				// 	for (int j = 0; j < nodesInSubnet.size(); j++) {
				// 		expr += X[i][j];
				// 	}
				// 	model.addConstr(I[i] <= 1);
				// 	model.addConstr(I[i] <= expr);
				// 	model.addConstr(nodesInSubnet.size() * I[i] >= expr);
				// }

				// Constraint 12
				patientsInCore = cores[subnetIdx].getSamplesPtr();
				for (int sampleIdx = 0; sampleIdx < patients.getSize(); sampleIdx++) {
					if (patientsInCore.getBit(sampleIdx) == 1) {
						model.addConstr(PTH[sourceIdx][sampleIdx] == 1);
					}
					else {
						model.addConstr(PTH[sourceIdx][sampleIdx] == 0);
					}
				}

				// Constraint 13
				for (int sampleIdx = 0; sampleIdx < patients.getSize(); sampleIdx++) {
					for (int j = 0; j < nodesInSubnet.size(); j++) {
						if (j != sourceIdx) {
							GRBLinExpr expr = 0;
							for (int i = 0; i < nodesInSubnet.size(); i++) {
								model.addConstr(PTHA[i][j][sampleIdx] <= PTH[i][sampleIdx]);
								model.addConstr(PTHA[i][j][sampleIdx] <= X[i][j]);
								model.addConstr(PTHA[i][j][sampleIdx] >= PTH[i][sampleIdx] + X[i][j] - 1);
								int helloAreYouThere = SubnetworkPatientsWithEdge[i][j]->getBit(sampleIdx);
								expr += helloAreYouThere * PTHA[i][j][sampleIdx];
							}
							model.addConstr(PTH[j][sampleIdx] == expr);
						}
					}
				}

				// Constraint 14
				for (int j = 0; j < nodesInSubnet.size(); j++) {
					if (j != sourceIdx) {
						GRBLinExpr expr = 0;
						for (int sampleIdx = 0; sampleIdx < patients.getSize(); sampleIdx++) {
							expr += PTH[j][sampleIdx];
						}
						// model.addConstr(expr >= minDepthRate);
						model.addConstr(expr >= double(ownPatients[j] -> getSize()) / 3.0);
					}
				}

				fprintf(stderr, "Finished constructing the model.\n");
				fprintf(stderr, "Solving ILP.\n");

				model.write("model.lp");
				model.optimize();

				char filename[1000];

				fprintf(foutTree, "***** %-30s *****\n", "Edge weights W_ij");
				for (int i = 0; i < nodesInSubnet.size(); i++) {
					for (int j = 0; j < nodesInSubnet.size(); j++) {
						fprintf(foutTree, "%.2lf ", W[i][j]);
					}
					fprintf(foutTree, "\n");
				}
				fprintf(foutTree, "\n");
				fprintf(foutTree, "adjM <- t(matrix(c(");
				for (int i = 0; i < nodesInSubnet.size(); i++) {
					for (int j = 0; j < nodesInSubnet.size(); j++) {
						if (i > 0 || j > 0) {
							fprintf(foutTree, ",");
						}
						fprintf(foutTree, "%.2lf", W[i][j]);
					}
				}
				fprintf(foutTree, "), nrow = %d, ncol = %d))\n\n", nodesInSubnet.size(), nodesInSubnet.size());

				fprintf(foutTree, "***** %-30s *****\n", "Thresholded Edge weights W_ij");
				for (int i = 0; i < nodesInSubnet.size(); i++) {
					for (int j = 0; j < nodesInSubnet.size(); j++) {
						// fprintf(foutTree, "%.2lf ", W[i][j]);
						if (W[i][j] >= edgeThreshold) {
							fprintf(foutTree, "%.2lf ", W[i][j]);
						}
						else {
							fprintf(foutTree, "0 ");
						}
					}
					fprintf(foutTree, "\n");
				}
				fprintf(foutTree, "\n");
				fprintf(foutTree, "adjM <- t(matrix(c(");
				for (int i = 0; i < nodesInSubnet.size(); i++) {
					for (int j = 0; j < nodesInSubnet.size(); j++) {
						if (i > 0 || j > 0) {
							fprintf(foutTree, ",");
						}
						// fprintf(foutTree, "%.2lf", W[i][j]);
						if (W[i][j] >= edgeThreshold) {
							fprintf(foutTree, "%.2lf", W[i][j]);
						}
						else {
							fprintf(foutTree, "0");
						}
					}
				}
				fprintf(foutTree, "), nrow = %d, ncol = %d))\n\n", nodesInSubnet.size(), nodesInSubnet.size());

				// fprintf(foutTree, "***** %-30s *****\n", "Patient agreement P_ij");
				// for (int i = 0; i < nodesInSubnet.size(); i++) {
				// 	for (int j = 0; j < nodesInSubnet.size(); j++) {
				// 		fprintf(foutTree, "%.2lf ", Agreed[i][j]);
				// 	}
				// 	fprintf(foutTree, "\n");
				// }
				// fprintf(foutTree, "\n");
				// fprintf(foutTree, "adjM <- t(matrix( c(");
				// for (int i = 0; i < nodesInSubnet.size(); i++) {
				// 	for (int j = 0; j < nodesInSubnet.size(); j++) {
				// 		if (i > 0 || j > 0) {
				// 			fprintf(foutTree, ", ");
				// 		}
				// 		fprintf(foutTree, "%.2lf", Agreed[i][j]);
				// 	}
				// }
				// fprintf(foutTree, "), nrow = %d, ncol = %d))\n\n", nodesInSubnet.size(), nodesInSubnet.size());

				fprintf(foutTree, "***** %-30s *****\n", "Precedence count adj_ij");
				for (int i = 0; i < nodesInSubnet.size(); i++) {
					for (int j = 0; j < nodesInSubnet.size(); j++) {
						fprintf(foutTree, "%d ", (int)adj[i][j]);
					}
					fprintf(foutTree, "\n");
				}
				fprintf(foutTree, "\n");
				fprintf(foutTree, "adjM <- t(matrix(c(");
				for (int i = 0; i < nodesInSubnet.size(); i++) {
					for (int j = 0; j < nodesInSubnet.size(); j++) {
						if (i > 0 || j > 0) {
							fprintf(foutTree, ",");
						}
						fprintf(foutTree, "%d", (int)adj[i][j]);
					}
				}
				fprintf(foutTree, "), nrow = %d, ncol = %d))\n\n", nodesInSubnet.size(), nodesInSubnet.size());

				fprintf(foutTree, "***** %-30s *****\n", "Tree edges X_ij - showing weight");
				for (int i = 0; i < nodesInSubnet.size(); i++) {
					for (int j = 0; j < nodesInSubnet.size(); j++) {
						// fprintf(foutTree, "%.1lf ", X[i][j].get(GRB_DoubleAttr_X));
						if (X[i][j].get(GRB_DoubleAttr_X) > 0.9) {
							fprintf(foutTree, "%.2lf ", W[i][j]);
						}
						else if (X[j][i].get(GRB_DoubleAttr_X) > 0.9) {
							fprintf(foutTree, "%.2lf ", -W[i][j]);
						}
						else {
							fprintf(foutTree, "0 ");
						}
					}
					fprintf(foutTree, "\n");
				}
				fprintf(foutTree, "\n");
				fprintf(foutTree, "adjM <- t(matrix(c(");
				for (int i = 0; i < nodesInSubnet.size(); i++) {
					for (int j = 0; j < nodesInSubnet.size(); j++) {
						if (i > 0 || j > 0) {
							fprintf(foutTree, ",");
						}
						// fprintf(foutTree, "%.1lf ", X[i][j].get(GRB_DoubleAttr_X));
						if (X[i][j].get(GRB_DoubleAttr_X) > 0.9) {
							fprintf(foutTree, "%.2lf", W[i][j]);
						}
						else if (X[j][i].get(GRB_DoubleAttr_X) > 0.9) {
							fprintf(foutTree, "%.2lf", -W[i][j]);
						}
						else {
							fprintf(foutTree, "0");
						}
					}
				}
				fprintf(foutTree, "), nrow = %d, ncol = %d))\n\n", nodesInSubnet.size(), nodesInSubnet.size());

				fprintf(foutTree, "***** %-30s *****\n", "Tree edges X_ij - showing precedence");
				for (int i = 0; i < nodesInSubnet.size(); i++) {
					for (int j = 0; j < nodesInSubnet.size(); j++) {
						// fprintf(foutTree, "%.1lf ", X[i][j].get(GRB_DoubleAttr_X));
						int numPaths = 0;
						for (int sampleIdx = 0; sampleIdx < patients.getSize(); sampleIdx++) {
							if (PTH[j][sampleIdx].get(GRB_DoubleAttr_X) > 0.9) {
								numPaths++;
							}
						}
						// int depth = D[j].get(GRB_DoubleAttr_X);
						// int numChildren = 0;
						// for (int k = 0; k < nodesInSubnet.size(); k++) {
						// 	if (X[j][k].get(GRB_DoubleAttr_X) > 0.9) {
						// 		numChildren++;
						// 	}
						// }
						if (X[i][j].get(GRB_DoubleAttr_X) > 0.9) {
							// fprintf(foutTree, "%d ", (int)adj[i][j]);
							fprintf(foutTree, "%d ", numPaths);
						}
						// else if (X[j][i].get(GRB_DoubleAttr_X) > 0.9) {
						// 	fprintf(foutTree, "%d ", -(int)adj[i][j]);
						// }
						else {
							fprintf(foutTree, "0 ");
						}
					}
					fprintf(foutTree, "\n");
				}
				fprintf(foutTree, "\n");
				fprintf(foutTree, "adjM <- t(matrix(c(");
				for (int i = 0; i < nodesInSubnet.size(); i++) {
					for (int j = 0; j < nodesInSubnet.size(); j++) {
						if (i > 0 || j > 0) {
							fprintf(foutTree, ",");
						}
						// fprintf(foutTree, "%.1lf ", X[i][j].get(GRB_DoubleAttr_X));
						int numPaths = 0;
						for (int sampleIdx = 0; sampleIdx < patients.getSize(); sampleIdx++) {
							if (PTH[j][sampleIdx].get(GRB_DoubleAttr_X) > 0.9) {
								numPaths++;
							}
						}
						// int depth = D[j].get(GRB_DoubleAttr_X);
						// int numChildren = 0;
						// for (int k = 0; k < nodesInSubnet.size(); k++) {
						// 	if (X[j][k].get(GRB_DoubleAttr_X) > 0.9) {
						// 		numChildren++;
						// 	}
						// }
						if (X[i][j].get(GRB_DoubleAttr_X) > 0.9) {
							// fprintf(foutTree, "%d", (int)adj[i][j]);
							fprintf(foutTree, "%d", numPaths);
						}
						// else if (X[j][i].get(GRB_DoubleAttr_X) > 0.9) {
						// 	fprintf(foutTree, "%d", -(int)adj[i][j]);
						// }
						else {
							fprintf(foutTree, "0");
						}
					}
				}
				fprintf(foutTree, "), nrow = %d, ncol = %d))\n\n", nodesInSubnet.size(), nodesInSubnet.size());

				// fprintf(foutTree, "***** %-30s *****\n", "Ancestor edges A_ij");
				// for (int i = 0; i < nodesInSubnet.size(); i++) {
				// 	for (int j = 0; j < nodesInSubnet.size(); j++) {
				// 		// fprintf(foutTree, "%.1lf ", X[i][j].get(GRB_DoubleAttr_X));
				// 		if (A[i][j].get(GRB_DoubleAttr_X) > 0.9) {
				// 			fprintf(foutTree, "%.2lf ", W[i][j]);
				// 		}
				// 		else {
				// 			fprintf(foutTree, "0 ");
				// 		}
				// 	}
				// 	fprintf(foutTree, "\n");
				// }
				// fprintf(foutTree, "\n");
				// fprintf(foutTree, "adjM <- t(matrix( c(");
				// for (int i = 0; i < nodesInSubnet.size(); i++) {
				// 	for (int j = 0; j < nodesInSubnet.size(); j++) {
				// 		if (i > 0 || j > 0) {
				// 			fprintf(foutTree, ", ");
				// 		}
				// 		// fprintf(foutTree, "%.1lf ", X[i][j].get(GRB_DoubleAttr_X));
				// 		if (A[i][j].get(GRB_DoubleAttr_X) > 0.9) {
				// 			fprintf(foutTree, "%.2lf ", W[i][j]);
				// 		}
				// 		else {
				// 			fprintf(foutTree, "0 ");
				// 		}
				// 	}
				// }
				// fprintf(foutTree, "), nrow = %d, ncol = %d))\n\n", nodesInSubnet.size(), nodesInSubnet.size());

				// fprintf(foutTree, "***** %30s *****\n", "Edge flow values");
				// for (int i = 0; i < nodesInSubnet.size(); i++) {
				// 	for (int j = 0; j < nodesInSubnet.size(); j++) {
				// 		fprintf(foutTree, "%.1lf ", F[i][j].get(GRB_DoubleAttr_X));
				// 	}
				// 	fprintf(foutTree, "\n");
				// }
				// fprintf(foutTree, "\n");
				// fprintf(foutTree, "adjM <- t(matrix( c(");
				// for (int i = 0; i < nodesInSubnet.size(); i++) {
				// 	for (int j = 0; j < nodesInSubnet.size(); j++) {
				// 		if (i > 0 || j > 0) {
				// 			fprintf(foutTree, ", ");
				// 		}
				// 		fprintf(foutTree, "%.1lf ", F[i][j].get(GRB_DoubleAttr_X));
				// 	}
				// }
				// fprintf(foutTree, "), nrow = %d, ncol = %d))\n\n", nodesInSubnet.size(), nodesInSubnet.size());

				// fprintf(foutTree, "***** %30s *****\n", "PD values");
				// for (int i = 0; i < nodesInSubnet.size(); i++) {
				// 	for (int j = 0; j < nodesInSubnet.size(); j++) {
				// 		fprintf(foutTree, "%.1lf ", PD[i][j].get(GRB_DoubleAttr_X));
				// 	}
				// 	fprintf(foutTree, "\n");
				// }
				// fprintf(foutTree, "adjM <- t(matrix( c(");
				// for (int i = 0; i < nodesInSubnet.size(); i++) {
				// 	for (int j = 0; j < nodesInSubnet.size(); j++) {
				// 		if (i > 0 || j > 0) {
				// 			fprintf(foutTree, ", ");
				// 		}
				// 		fprintf(foutTree, "%.1lf ", PD[i][j].get(GRB_DoubleAttr_X));
				// 	}
				// }
				// fprintf(foutTree, "), nrow = %d, ncol = %d))\n\n", nodesInSubnet.size(), nodesInSubnet.size());

				// fprintf(foutTree, "***** %-30s *****\n", "Inner node flag");
				// for (int i = 0; i < nodesInSubnet.size(); i++) {
				// 	int val = I[i].get(GRB_DoubleAttr_X);
				// 	fprintf(foutTree, "%d [%s : %d]\n", i + 1, nodes.names[nodeIdxInSubnet[i]].c_str(), val);
				// }
				// fprintf(foutTree, "\n");

				fprintf(foutTree, "***** %-30s *****\n", "Node depth");
				for (int i = 0; i < nodesInSubnet.size(); i++) {
					int depth = D[i].get(GRB_DoubleAttr_X);
					fprintf(foutTree, "%d [%s : %d]\n", i + 1, nodes.names[nodeIdxInSubnet[i]].c_str(), depth);
				}
				fprintf(foutTree, "\n");

				fprintf(foutTree, "qgraph(adjM, edge.labels = TRUE, labels = c(");
				for (int i = 0; i < nodesInSubnet.size(); i++) {
					if (i > 0) {
						fprintf(foutTree, ", ");
					}
					// int numPaths = 0;
					// for (int sampleIdx = 0; sampleIdx < patients.getSize(); sampleIdx++) {
					// 	if (PTH[i][sampleIdx].get(GRB_DoubleAttr_X) > 0.9) {
					// 		numPaths++;
					// 	}
					// }
					fprintf(foutTree, "\"%d %s\"", (int)adj[sourceIdx][i], nodes.names[nodeIdxInSubnet[i]].c_str());
					// fprintf(foutTree, "\"%d %s\"", numPaths, nodes.names[nodeIdxInSubnet[i]].c_str());
				}
				fprintf(foutTree, "))\n\n");

				pair<int, int> * nodeDepths = new pair<int, int> [nodesInSubnet.size()];
				int * nodeStack = new int[nodesInSubnet.size()];
				for (int i = 0; i < nodesInSubnet.size(); i++) {
					int depth = D[i].get(GRB_DoubleAttr_X);
					nodeDepths[i] = make_pair(i, depth);
				}
				sort(nodeDepths, nodeDepths + nodesInSubnet.size(), cmpPair);
				for (int i = nodesInSubnet.size() - 1; i >= 0; i--) {
					int nodeIdx = nodeDepths[i].first;
					int depth = nodeDepths[i].second;
					int stackSize = 0;
					while (depth--) {
						nodeStack[stackSize++] = nodeIdx;
						int parentIdx = -1;
						for (int j = nodesInSubnet.size() - 1; j >= 0; j--) {
							if (X[j][nodeIdx].get(GRB_DoubleAttr_X) > 0.9) {
								parentIdx = j;
								break;
							}
						}
						nodeIdx = parentIdx;
					}
					fprintf(foutTree, "%s", nodes.names[nodeIdxInSubnet[sourceIdx]].c_str());
					while (stackSize > 0) {
						int nodeIdx = nodeStack[--stackSize];
						int numPaths = 0;
						for (int sampleIdx = 0; sampleIdx < patients.getSize(); sampleIdx++) {
							if (PTH[nodeIdx][sampleIdx].get(GRB_DoubleAttr_X) > 0.9) {
								numPaths++;
							}
						}
						fprintf(foutTree, " --[%d]--> %s (%d)", numPaths, nodes.names[nodeIdxInSubnet[nodeIdx]].c_str(), (int)adj[sourceIdx][nodeIdx]);
					}
					fprintf(foutTree, "\n");
				}
				fclose(foutTree);

				delete nodeDepths;
				delete nodeStack;
				/**/
			}
			catch (GRBException e){
				cout << "Error code = " << e.getErrorCode() << endl;
				cout << e.getMessage() << endl;
			}
			catch (...) {
				cout << "Exception during optimization" << endl;
			}

			// Freeing the memory
			for (int i = 0; i < nodesInSubnet.size(); i++) {
				if (F && F[i])
					delete[] F[i];
				if (X && X[i])
					delete[] X[i];
				if (PD && PD[i])
					delete[] PD[i];
			}
			delete F;
			delete X;
			delete PD;
			delete [] D;
			// delete ub;
			delete env;
			for (auto it : adjIdx) {
				int idxInMatrix = it.second;
				delete adj[idxInMatrix];
				delete W[idxInMatrix];
				delete Agreed[idxInMatrix];
			}
			delete W;
			delete Agreed;
			delete adj;
			for (int nodeIdx = 0; nodeIdx < nodesInSubnet.size(); nodeIdx++) {
				for (int neighbourIdx = 0; neighbourIdx < nodesInSubnet.size(); neighbourIdx++) {
					delete SubnetworkPatientsWithEdge[nodeIdx][neighbourIdx];
				}
				delete SubnetworkPatientsWithEdge[nodeIdx];
				delete ownPatients[nodeIdx];
			}
			delete SubnetworkPatientsWithEdge;
			delete ownPatients;
		}
	}

	// Deallocating memory
	for (int altIdx = 0; altIdx < alterations.getSize(); altIdx++) {
		for (int nodeIdx = 0; nodeIdx < G.V; nodeIdx++) {
			delete overlap[altIdx][nodeIdx];
		}
		delete overlap[altIdx];
	}
	delete overlap;

	return maxSize;
}

void calculatePValues_permutationTest(Entry & patients, Entry & nodes, Entry & alterations, Bitmask *** geneSampleAlterations, Graph & G, unordered_map<string, Bitmask> & edgePatients, string outFolder, const int minDepth, const double recurrenceRate, const int numIterations, const int maxFoundSubnetSize, bool permuteOnlyOrder = false, const char * nodeSetFilename = NULL) {
	// int numGreaterThanOrEqual = 0;
	int * numGreaterThanOrEqual = new int[maxFoundSubnetSize + 1];
	memset(numGreaterThanOrEqual, 0, sizeof(numGreaterThanOrEqual[0]) * (maxFoundSubnetSize + 1));
	char filename[1000];
	sprintf(filename, "%s/permutationSizeDistribution.txt", outFolder.c_str());
	FILE * foutDist = fopen(filename, "w");
	int ** numEvents = new int * [patients.getSize()];
	vector< pair<int, int> > * originalNodes = new vector< pair<int, int> > [patients.getSize()];
	vector<string> * originalEdges = new vector<string> [patients.getSize()];
	random_device rd;
	mt19937_64 * RNG;
	RNG = new std::mt19937_64(rd());
	for (int sampleIdx = 0; sampleIdx < patients.getSize(); sampleIdx++) {
		numEvents[sampleIdx] = new int [alterations.getSize()];
		memset(numEvents[sampleIdx], 0, sizeof(numEvents[sampleIdx][0]) * alterations.getSize());
	}
	Bitmask *** new_geneSampleAlterations = new Bitmask ** [alterations.getSize()];
	for (int altIdx = 0; altIdx < alterations.getSize(); altIdx++) {
		new_geneSampleAlterations[altIdx] = new Bitmask * [nodes.getSize()];
		for (int nodeIdx = 0; nodeIdx < nodes.getSize(); nodeIdx++) {
			new_geneSampleAlterations[altIdx][nodeIdx] = nullptr;
			if (geneSampleAlterations[altIdx][nodeIdx] != nullptr) {
				Bitmask tempMask(geneSampleAlterations[altIdx][nodeIdx]);
				while (tempMask.getSize()) {
					int sampleIdx = tempMask.extractLowestOrderSetBitIndex();
					numEvents[sampleIdx][altIdx]++;
					originalNodes[sampleIdx].push_back( make_pair(nodeIdx, altIdx) );
				}
			}
		}
	}
	for (auto it : edgePatients) {
		Bitmask tempMask(it.second);
		while (tempMask.getSize()) {
			int sampleIdx = tempMask.extractLowestOrderSetBitIndex();
			originalEdges[sampleIdx].push_back(it.first);
		}
	}

	// Generating random alteration permutations across patients
	for (int it_idx = 0, lastProg = 0; it_idx < numIterations; it_idx++) {
		int progress = 1000 * double(it_idx + 1) / double(numIterations);
		if (progress > lastProg) {
			fprintf(stderr, "\r\t%.1lf%%", double(progress)/10.0);
			lastProg = progress;
		}
		unordered_map<string, Bitmask> new_edgePatients;
		Graph new_G;
		unordered_set<pair<string, string>, stringPairHash> uniqueEdges;
		for (int altIdx = 0; altIdx < alterations.getSize(); altIdx++) {
			for (int nodeIdx = 0; nodeIdx < nodes.getSize(); nodeIdx++) {
				delete new_geneSampleAlterations[altIdx][nodeIdx];
				new_geneSampleAlterations[altIdx][nodeIdx] = nullptr;
			}
		}
		bool weGoAgane = false;
		for (int sampleIdx = 0; sampleIdx < patients.getSize(); sampleIdx++) {
			// (Step 1) Choosing random nodes/events for the current patient (preserving the original event type count)
			vector<int> eventNodes;
			vector<int> eventColors;
			unordered_set<string> CNAAlreadyPicked;
			// We want to leave the germline root in the same place. For now, add it first.
			eventNodes.push_back(nodes.indices["GL"]);
			eventColors.push_back(alterations.indices["-"]);


			for (int altIdx = 0; altIdx < alterations.getSize(); altIdx++) {
				if (alterations.names[altIdx] == "-") {
					continue;
				}
				vector<int> eventPool;
				vector<int> eventRecurrence;
				llu totalRecurrence = 0;
				for (int nodeIdx = 0; nodeIdx < nodes.getSize(); nodeIdx++) {
					if (geneSampleAlterations[altIdx][nodeIdx] != nullptr && geneSampleAlterations[altIdx][nodeIdx]->getSize() > 0 && !CNAAlreadyPicked.count(nodes.names[nodeIdx])) {
						eventPool.push_back(nodeIdx);
						eventRecurrence.push_back(geneSampleAlterations[altIdx][nodeIdx]->getSize());
						totalRecurrence += geneSampleAlterations[altIdx][nodeIdx]->getSize();
					}
				}
				int eventPoolSize = eventPool.size();
				int numColors = numEvents[sampleIdx][altIdx];
				// fprintf(stderr, "Event pool: %d\t\tNumColors: %d\n", eventPoolSize, numColors);
				
				while (numColors > 0) {
					// Choose an event
					if (eventPoolSize == 0) {
						fprintf(stderr, "\n< Error > We ran out of events... Awkward... Trying again!\n");
						weGoAgane = true;
						// exit(0);
						break;
					}
					std::uniform_int_distribution<int> uni(0, totalRecurrence);
					int r = uni(*RNG);
					int randIdx = 0;
					while (r > eventRecurrence[randIdx]) {
						r-= eventRecurrence[randIdx];
						randIdx++;
					}
					int nodeIdx = eventPool[randIdx];
						// for (int i = 0; i < eventPoolSize; i++) {
						// 	fprintf(stderr, "%s\t", nodes.names[eventPool[i]].c_str());
						// }
						// fprintf(stderr, "\n> ");
					// Remove the event so we can't draw it again
					totalRecurrence -= eventRecurrence[randIdx];
					swap(eventPool[randIdx], eventPool[eventPoolSize - 1]);
					swap(eventRecurrence[randIdx], eventRecurrence[eventPoolSize - 1]);
					eventPoolSize--;
					// If it was CNA, we have checks to perform
					bool thisIsWrong = false;
					if (alterations.names[altIdx] == "CNLOSS" || alterations.names[altIdx] == "CNGAIN") {
						// fprintf(stderr, "%s\n", nodes.names[nodeIdx].c_str());
						string chrLabel = nodes.names[nodeIdx].substr(5);
						char chrArm = chrLabel[chrLabel.size() - 1];
						int chrNum;
						sscanf(chrLabel.c_str(), "%d", &chrNum);
						char fullLabel[1000];
						// If there was a whole-chromosome event on the same chromsome in this patient, then we can't assign the current event
						sprintf(fullLabel, "gain_%d", chrNum);
						thisIsWrong |= nodes.indices.count(fullLabel) && CNAAlreadyPicked.count(fullLabel);
						sprintf(fullLabel, "loss_%d", chrNum);
						thisIsWrong |= nodes.indices.count(fullLabel) && CNAAlreadyPicked.count(fullLabel);
						// If it's an arm event, we can't already have a loss followed by a gain on the same arm (cannot be controlled in this code segment)
						if (chrArm == 'p' || chrArm == 'q') {
							// sprintf(fullLabel, "gain_%d%c", chrNum, chrArm);
							// thisIsWrong |= nodes.indices.count(fullLabel) && CNAAlreadyPicked.count(fullLabel);
							// sprintf(fullLabel, "loss_%d%c", chrNum, chrArm);
							// thisIsWrong |= nodes.indices.count(fullLabel) && CNAAlreadyPicked.count(fullLabel);
						}
						// If it's a whole chromosome event, we can't have either arm event in the same patient
						else {
							sprintf(fullLabel, "gain_%dp", chrNum);
							thisIsWrong |= nodes.indices.count(fullLabel) && CNAAlreadyPicked.count(fullLabel);
							sprintf(fullLabel, "gain_%dq", chrNum);
							thisIsWrong |= nodes.indices.count(fullLabel) && CNAAlreadyPicked.count(fullLabel);
							sprintf(fullLabel, "loss_%dp", chrNum);
							thisIsWrong |= nodes.indices.count(fullLabel) && CNAAlreadyPicked.count(fullLabel);
							sprintf(fullLabel, "loss_%dq", chrNum);
							thisIsWrong |= nodes.indices.count(fullLabel) && CNAAlreadyPicked.count(fullLabel);
						}
					}
					// fprintf(stderr, "\n123\n");
					// If everything is fine, add it to the event list
					if (!thisIsWrong) {
						eventNodes.push_back(nodeIdx);
						eventColors.push_back(altIdx);
						if (alterations.names[altIdx] == "CNLOSS" || alterations.names[altIdx] == "CNGAIN") {
							CNAAlreadyPicked.insert( nodes.names[nodeIdx] );
						}
						numColors--;
						// fprintf(stderr, "%s", nodes.names[nodeIdx].c_str());
					}
					// fprintf(stderr, "\n");
				}
				if (weGoAgane) {
					break;
				}
			}
			if (weGoAgane) {
				break;
			}

			// (Step 2) Generating a random order of the original nodes.
			vector< pair<int, int> > permutedNodes(originalNodes[sampleIdx]);
			for (int i = 0; i < permutedNodes.size() - 1; i++) {
				std::uniform_int_distribution<int> uni(i, permutedNodes.size() - 1);
				int randIdx = uni(*RNG);
				if (i != randIdx) {
					swap(permutedNodes[i], permutedNodes[randIdx]);
				}
			}
			// Put the germline root as first so that it matches GL in the randomly chosen node set
			int GLIdx = 0;
			while (GLIdx < permutedNodes.size() && permutedNodes[GLIdx].first != nodes.indices["GL"]) {
				GLIdx++;
			}
			if (GLIdx > 0) {
				swap(permutedNodes[0], permutedNodes[GLIdx]);
			}

			// (Step 3) Renaming edges.
			unordered_map<string, string> nodeNameMap;
			for (int i = 0; i < permutedNodes.size() - 1; i++) {
				string node1 = nodes.names[eventNodes[i]];
				string alt1 = alterations.names[eventColors[i]];
				string node2 = nodes.names[permutedNodes[i].first];
				string alt2 = alterations.names[permutedNodes[i].second];
				string originalEventName = node2 + " " + alt2;
				string newEventName = node1 + " " + alt1;
				nodeNameMap[originalEventName] = newEventName;
			}
			for (string & originalEdgeName : originalEdges[sampleIdx]) {
				char node1[1000];
				char alt1[1000];
				char node2[1000];
				char alt2[1000];
				sscanf(originalEdgeName.c_str(), "%s %s %s %s", node1, alt1, node2, alt2);
				string event1 = nodeNameMap[string(node1) + " " + string(alt1)];
				string event2 = nodeNameMap[string(node2) + " " + string(alt2)];
				string newEdgeName = event1 + " " + event2;
				new_edgePatients.emplace(newEdgeName, patients.getSize());
				new_edgePatients.find(newEdgeName)->second.setBit(sampleIdx, 1);
				sscanf(event1.c_str(), "%s", node1);
				sscanf(event2.c_str(), "%s", node2);
				pair<string, string> e = make_pair(string(node1), string(node2));
				uniqueEdges.insert(e);
			}

			// (Step 4) Filling out new patient alteration bitmasks
			for (int i = 0; i < eventNodes.size(); i++) {
				int nodeIdx = eventNodes[i];
				int altIdx = eventColors[i];
				if (new_geneSampleAlterations[altIdx][nodeIdx] == nullptr) {
					new_geneSampleAlterations[altIdx][nodeIdx] = new Bitmask(patients.getSize());
				}
				new_geneSampleAlterations[altIdx][nodeIdx]->setBit(sampleIdx, 1);
			}
		}
		if (weGoAgane) {
			it_idx--;
			continue;
		}

		// (Step 5) Generating new graph new_G from new_edgePatients
		new_G.E = uniqueEdges.size();
		new_G.V = nodes.indices.size();
		new_G.NSize = new int [G.V];
		new_G.N = new int * [G.V];
		memset(new_G.NSize, 0, sizeof(new_G.NSize[0]) * new_G.V);
		for (auto e: uniqueEdges) {
			int idx1 = nodes.indices[e.first];
			new_G.NSize[idx1]++;
		}
		for (int i = 0; i < new_G.V; i++) {
			new_G.N[i] = new int[ new_G.NSize[i] ];
		}
		int * tempNSize = new int [new_G.V];
		memset(tempNSize, 0, sizeof(tempNSize[0]) * new_G.V );
		for (auto e: uniqueEdges) {
			int idx1 = nodes.indices[e.first];
			int idx2 = nodes.indices[e.second];
			new_G.N[ idx1 ][ tempNSize[idx1] ] = idx2;
			tempNSize[idx1]++;
		}
		delete [] tempNSize;

		// Computing the subnetworks
		llu numSubnetworks = 0;
		SubnetworkEntry * coreSubnetworks = nullptr;
		// constructSingleNodeSubnetworks(minDepth, coreSubnetworks, numSubnetworks, patients, nodes, alterations, new_geneSampleAlterations, false);
		if (nodeSetFilename != NULL) {
			constructSubnetworkCoresFromSet(nodeSetFilename, coreSubnetworks, numSubnetworks, patients, nodes, alterations, new_geneSampleAlterations, false);
		}
		else {
			constructSingleNodeSubnetworks(minDepth, coreSubnetworks, numSubnetworks, patients, nodes, alterations, new_geneSampleAlterations, false);
		}
		// fprintf(stderr, "\n\n> ");
		// for (int i = 0; i < numSubnetworks; i++) {
		// 	const map<int, llu> & nodesInSubnet = coreSubnetworks[i].getNodesRef();
		// 	for (auto it : nodesInSubnet) {
		// 		int nodeIdx = it.first;
		// 		llu altMask = it.second;
		// 		fprintf(stderr, "%s\t", nodes.names[nodeIdx].c_str());
		// 	}
		// }
		// fprintf(stderr, "\n");
		SubnetworkEntry * candidateSubnetworks = new SubnetworkEntry [numSubnetworks];
		constructSubnetworks(coreSubnetworks, numSubnetworks, recurrenceRate, candidateSubnetworks, patients, nodes, alterations, new_geneSampleAlterations, new_G, new_edgePatients, false);
		// constructSubnetworks(coreSubnetworks, numSubnetworks, minDepth, candidateSubnetworks, patients, nodes, alterations, new_geneSampleAlterations, new_G, new_edgePatients, false);
		// fprintf(stderr, "\n\n> ");
		// for (int i = 0; i < numSubnetworks; i++) {
		// 	{
		// 		const map<int, llu> & nodesInSubnet = coreSubnetworks[i].getNodesRef();
		// 		for (auto it : nodesInSubnet) {
		// 			int nodeIdx = it.first;
		// 			llu altMask = it.second;
		// 			fprintf(stderr, "%s\t", nodes.names[nodeIdx].c_str());
		// 		}
		// 	}
		// 	const map<int, llu> & nodesInSubnet = candidateSubnetworks[i].getNodesRef();
		// 	fprintf(stderr, " (");
		// 	for (auto it : nodesInSubnet) {
		// 		int nodeIdx = it.first;
		// 		llu altMask = it.second;
		// 		fprintf(stderr, "%s\t", nodes.names[nodeIdx].c_str());
		// 	}
		// 	fprintf(stderr, " )");
		// }
		// fprintf(stderr, "\n");
		// Finding the maximum size
		int maxSize = 0;
		if (numSubnetworks > 0) {
			int maxIdx = 0;
			int maxCnt = 1;
			for (int i = 1; i < numSubnetworks; i++) {
				if (candidateSubnetworks[i].getSize() > candidateSubnetworks[maxIdx].getSize()) {
					maxIdx = i;
					maxCnt = 1;
				}
				else if (candidateSubnetworks[i].getSize() == candidateSubnetworks[maxIdx].getSize()) {
					maxCnt++;
				}
			}
			maxSize = candidateSubnetworks[maxIdx].getSize();
			// for (int i = 0; i < numSubnetworks; i++) {
			// 	if (candidateSubnetworks[i].getSize() == maxSize) {
			// 		const map<int, llu> & nodesInSubnet = candidateSubnetworks[i].getNodesRef();
			// 		fprintf(stderr, "\n> ");
			// 		for (auto it : nodesInSubnet) {
			// 			int nodeIdx = it.first;
			// 			llu altMask = it.second;
			// 			fprintf(stderr, "%s", nodes.names[nodeIdx].c_str());
			// 			Bitmask alterationMask(64);
			// 			alterationMask.copylluBitmask(altMask);	// Go through alterations of the node
			// 			while (alterationMask.getSize()) {
			// 				int altIdx = alterationMask.extractLowestOrderSetBitIndex();
			// 				fprintf(stderr, " (%s)\t", alterations.names[altIdx].c_str());
			// 			}
			// 		}
			// 		fprintf(stderr, "\n");
			// 	}
			// }
		}
		// fprintf(stderr, "Maximum subnetwork size is %d.\n", maxSize);
		delete [] coreSubnetworks;
		delete [] candidateSubnetworks;

		// if (maxSize >= maxFoundSubnetSize) {
		// 	numGreaterThanOrEqual++;
		// }
		if (maxSize > maxFoundSubnetSize) {
			maxSize = maxFoundSubnetSize;
		}
		for (int i = maxSize; i > 0; i--) {
			numGreaterThanOrEqual[i]++;
		}
		fprintf(foutDist, "%d\n", maxSize);
	}

	for (int altIdx = 0; altIdx < alterations.getSize(); altIdx++) {
		for (int nodeIdx = 0; nodeIdx < nodes.getSize(); nodeIdx++) {
			delete new_geneSampleAlterations[altIdx][nodeIdx];
		}
		delete new_geneSampleAlterations[altIdx];
	}
	delete new_geneSampleAlterations;
	for (int sampleIdx = 0; sampleIdx < patients.getSize(); sampleIdx++) {
		delete numEvents[sampleIdx];
	}
	delete numEvents;
	delete [] originalNodes;
	delete [] originalEdges;
	delete RNG;

	fclose(foutDist);
	// fprintf(stderr, "\nRandom event assignment in each patient results in the subnetwork of size at least %d in %d / %d experiments.\n", maxFoundSubnetSize, numGreaterThanOrEqual, numIterations);
	for (int i = 2; i <= maxFoundSubnetSize; i++) {
		fprintf(stderr, "\nRandom event assignment in each patient results in the subnetwork of size at least %d in %d / %d experiments.\n", i, numGreaterThanOrEqual[i], numIterations);
	}
}

int main(int argc, char * argv[]) {
	printHeader("CONETT");
	// INPUT CHECK
	if (argc <= 1) {
		fprintf(stderr, "./cd-RECAP_phylogeny -p [input file with phylogenies] -g [gene sets to be used as cores] -t [minimum seed recurrence] -e [subgraph tumor conservation rate] -i [permutation p-value iterations] -a [edge weight cutoff] -f [outputFolder] \n\n");
		return 0;
	}
	char consoleFlags[] = {'p', 'e', 't', 'f', 'g', 'i', 'a', 0};
	bool optional[200] = {};
	optional['g'] = true;
	optional['f'] = true;
	optional['i'] = true;
	optional['a'] = true;
	unordered_map<char, string> consoleParameters;
	for (int i = 1; i < argc; i++) {
		if ( argv[i][0] == '-' && argv[i][1] && i + 1 < argc && argv[i + 1][0] != '-' ) {
			consoleParameters[ argv[i][1] ] = string( argv[i + 1] );
			i++;
		}
	}
	for (char * ptrFlag = consoleFlags; *ptrFlag; ptrFlag++) {
		if ( !optional[*ptrFlag] && !consoleParameters.count(*ptrFlag) ) {
			fprintf(stderr, "\n< Error > Missing value for parameter '%c'. Exiting program.\n", *ptrFlag);
			exit(0);
		}
	}
	if (!consoleParameters.count('f')) {
		consoleParameters['f'] = "testRun";
		fprintf(stderr, "Value for parameter 'f' not specified. Writing all output to \"testRun\" folder.\n");
	}
	if (!consoleParameters.count('i')) {
		consoleParameters['i'] = "0";
	}
	if (!consoleParameters.count('a')) {
		consoleParameters['a'] = "0.85";
	}
	// READING INPUT PARAMETERS
	int minDepth;
	double edgeThreshold;
	double recurrenceRate;
	int pValPermutationIterations;
	char folderName[1000] = {};
	sscanf(consoleParameters['t'].c_str(), "%d", &minDepth);
	sscanf(consoleParameters['i'].c_str(), "%d", &pValPermutationIterations);
	sscanf(consoleParameters['e'].c_str(), "%lf", &recurrenceRate);
	sscanf(consoleParameters['a'].c_str(), "%lf", &edgeThreshold);
	sscanf(consoleParameters['f'].c_str(), "%s", folderName);
	// CREATING DIRECTORY STRUCTURE FOR THE OUTPUT
	char command[1000];
	char fullFolder[1000];
	sprintf(fullFolder, "./%s_t%d_e%2lf", ("output/" + string(folderName)).c_str(), minDepth, recurrenceRate);
	string outFolder = string(fullFolder);
	fprintf(stderr, "Using output folder %s\n", outFolder.c_str());
	sprintf(command, "rm -f -r %s", outFolder.c_str());
	system(command);
	sprintf(command, "mkdir -p %s", outFolder.c_str());
	system(command);
	// CREATING DATA STRUCTURES
	Entry patients, nodes, alterations;
	Graph G;
	Bitmask *** geneSampleAlterations;
	unordered_map<string, Bitmask> edgePatients;
	// unordered_map<int, llu> * geneAlterationBitmasks;	// geneAlterationBitmasks[ i ][ j ] = c means that "gene i has colour c in patient j". The colours are bitmasks (so supporting max 64 different alteration types).
	// unordered_map<string, string> geneChrLoc;
	// READING INPUT FILES
	printHeader("Processing Input");
	readAlterationPhylogeny(consoleParameters['p'].c_str(), patients, nodes, alterations, geneSampleAlterations, G, edgePatients, outFolder.c_str());
	// CONSTRUCTING CORE SUBNETWORKS
	llu numCoreSubnetworks = 0;
	SubnetworkEntry * coreSubnetworks = nullptr;
	if (consoleParameters.count('g')) {
		printHeader("Constructing core subnetworks from gene sets");
		constructSubnetworkCoresFromSet(consoleParameters['g'].c_str(), coreSubnetworks, numCoreSubnetworks, patients, nodes, alterations, geneSampleAlterations);
	}
	else {
		printHeader("Constructing core subnetworks of single genes");
		constructSingleNodeSubnetworks(minDepth, coreSubnetworks, numCoreSubnetworks, patients, nodes, alterations, geneSampleAlterations);
	}
	if (numCoreSubnetworks == 0)
		return 0;
	// EXTENDING CORE SUBNETWORKS
	printHeader("Extending core subnetworks");
	llu numSubnetworks = numCoreSubnetworks;
	SubnetworkEntry * candidateSubnetworks = new SubnetworkEntry [numSubnetworks];
	constructSubnetworks(coreSubnetworks, numCoreSubnetworks, recurrenceRate, candidateSubnetworks, patients, nodes, alterations, geneSampleAlterations, G, edgePatients);
	// constructSubnetworks(coreSubnetworks, numCoreSubnetworks, minDepth, candidateSubnetworks, patients, nodes, alterations, geneSampleAlterations, G, edgePatients);
	// for (int i = 0; i < numSubnetworks; i++) {
	// 	fprintf(stderr, "#%d size: %d\n", i, candidateSubnetworks[i].getSize());
	// 	const map<int, llu> & nodesInSubnet = candidateSubnetworks[i].getNodesRef();
	// 	for (auto it : nodesInSubnet) {
	// 		int nodeIdx = it.first;
	// 		llu altMask = it.second;
	// 		fprintf(stderr, "%s\t", nodes.names[nodeIdx].c_str());
	// 	}
	// 	fprintf(stderr, "\n");
	// }
	printHeader("Choosing the optimal subnetwork");
	int maxFoundSubnetSize = chooseOptimalSubnetwork(candidateSubnetworks, numSubnetworks, recurrenceRate, edgeThreshold, coreSubnetworks, outFolder, patients, nodes, alterations, geneSampleAlterations, G, edgePatients);
	// int maxFoundSubnetSize = chooseOptimalSubnetwork(candidateSubnetworks, numSubnetworks, minDepth, edgeThreshold, coreSubnetworks, outFolder, patients, nodes, alterations, geneSampleAlterations, G, edgePatients);
	if (pValPermutationIterations > 0) {
		printHeader("Calculating the p-value via random node event reassignment test");
		calculatePValues_permutationTest(patients, nodes, alterations, geneSampleAlterations, G, edgePatients, outFolder, minDepth, recurrenceRate, pValPermutationIterations, maxFoundSubnetSize, false, (consoleParameters.count('g') ? consoleParameters['g'].c_str() : NULL));
	}
	// int t[] = {5, 6, 7, 8, 9};
	// int numt = sizeof(t) / sizeof(t[0]);
	// for (int i = 0; i < numt; i++) {
	// // for (int i = 7; i < numt; i++) {
	// 	for (double consR = 0.8; consR >= 0.5; consR -= 0.1) {
	// 		fprintf(stderr, "t = %d\tconsR = %.1lf\n", t[i], consR);
	// 		calculatePValues_permutationTest(patients, nodes, alterations, geneSampleAlterations, G, edgePatients, outFolder, t[i], consR, pValPermutationIterations, 10, false);
	// 	}
	// }
	return 0;
}