#ifndef kmer_hpp
#define kmer_hpp
#include <fstream>
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "sequence_reverse.hpp"

using namespace std;

vector<string> build_kmer_index(const string & sequence, const int kmerLen);
string sequence_select(const string & sequence);
float get_jacccard(const vector<string> & minimizerVectorBase, const vector<string> & minimizerVectorQry);
float get_base_jacccard(const vector<string> & minimizerVectorBase, const vector<string> & minimizerVectorQry);
float get_qry_jacccard(const vector<string> & minimizerVectorBase, const vector<string> & minimizerVectorQry);


vector<string> build_kmer_index(const string & sequence, const int kmerLen) {
    int seqLen = sequence.size();
    int kmerNumber = seqLen - kmerLen + 1;
    
    int seqStart = 0;
    // Vector deduplication
    vector<string> kmerVector;

    for (int i = 0; i < kmerNumber; i++) {
        // Extract kmer
        string kmerSeq = sequence.substr(seqStart, kmerLen);

        // Select the smallest sequence between the reverse complement sequence
        kmerSeq = sequence_select(kmerSeq);

        // Add to kmerVector
        kmerVector.push_back(kmerSeq);

        seqStart += 1;
    }

    return kmerVector;
}


// Reverse complement of a sequence and return the minimum value
string sequence_select(const string & sequence) {
    // Sequence reverse complement
    string sequenceRev = sequence_reverse(sequence);;
    
    // If it contains special characters, the original string is returned.
    if (sequenceRev.empty()) {
        return sequence;
    }
    
    // Otherwise, return the smallest string.
    return min(sequence, sequenceRev);
}


// Get the Jaccard similarity between sequences
float get_qry_jacccard(const vector<string> & minimizerVectorBase, const vector<string> & minimizerVectorQry) {
    float jacccard;
    int matchMinimizerNumber = 0;
    for (auto it : minimizerVectorQry) {
        if (find(minimizerVectorBase.begin(), minimizerVectorBase.end(), it) != minimizerVectorBase.end()) {
            matchMinimizerNumber += 1;
        }
    }
    
    jacccard = matchMinimizerNumber/(float)minimizerVectorQry.size();
    return jacccard;
}

float get_base_jacccard(const vector<string> & minimizerVectorBase, const vector<string> & minimizerVectorQry) {
    float jacccard;
    int matchMinimizerNumber = 0;
    for (auto it : minimizerVectorQry) {
        if (find(minimizerVectorBase.begin(), minimizerVectorBase.end(), it) != minimizerVectorBase.end()) {
            matchMinimizerNumber += 1;
        }
    }
    
    jacccard = matchMinimizerNumber/(float)minimizerVectorBase.size();
    return jacccard;
}

float get_jacccard(const vector<string> & minimizerVectorBase, const vector<string> & minimizerVectorQry) {
    float jacccard;
    int matchMinimizerNumber = 0;
    for (auto it : minimizerVectorQry) {
        if (find(minimizerVectorBase.begin(), minimizerVectorBase.end(), it) != minimizerVectorBase.end()) {
            matchMinimizerNumber += 1;
        }
    }
    
    jacccard = matchMinimizerNumber/(float)minimizerVectorBase.size();
    return jacccard;
}

#endif /* kmer_hpp */