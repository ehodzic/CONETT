#ifndef __ENTRY_H_INCLUDED__
#define __ENTRY_H_INCLUDED__
#include <unordered_map>

struct Entry {
	std::string * names;
	std::unordered_map<std::string, int> indices;

	Entry();

	int getSize() const { return indices.size(); }
};

#endif