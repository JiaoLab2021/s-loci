#include <fstream>
#include <string>
#include <algorithm>
#include <unordered_map>

using namespace std;

unordered_map<char, char> ntRevMap = {
	{'a', 't'}, {'t', 'a'}, {'g', 'c'}, {'c', 'g'}, {'n', 'n'}, 
	{'A', 'T'}, {'T', 'A'}, {'G', 'C'}, {'C', 'G'}, {'N', 'N'}
};

// Reverse complement
string sequence_reverse(const string & sequence) {
	// Reverse complement
	string sequenceRev = sequence;
	reverse(sequenceRev.begin(), sequenceRev.end());

	for (int i = 0; i < sequenceRev.size(); i++) {
		auto findIter = ntRevMap.find(sequenceRev[i]);
		if (findIter != ntRevMap.end()) {
			sequenceRev[i] = findIter->second;
		} else {
			sequenceRev = "";
			return sequenceRev;
		}
	}

	return sequenceRev;
}