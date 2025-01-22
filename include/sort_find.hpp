#ifndef sort_find_hpp
#define sort_find_hpp
#include <iostream>
#include <algorithm>
#include <vector>
#include <regex>
#include <limits.h>

using namespace std;

template<typename T>

// Array Indexing
int findPosArray(T ar[], int n, T element) {  // Find an element and return the position index, find(array, length, element)
	int i = 0;
	int index=-1;  // Original subscript, returns -1 if no element is found
	for (i = 0; i <n; i++) {
		if (element ==ar[i]) {
			index=i;  // Record element subscript
		}
	}
	return index;  // Return subscript
}

// Vector Indexing
int findPosVector(vector <int> input , int number) {
    vector<int>::iterator iter=std::find(input.begin(),input.end(),number);  // What is returned is an iterator pointer
    if(iter == input.end()) {
        return -1;
    } else {
        return std::distance(input.begin(),iter);
    }
}

// Reverse
bool Reverse(int a,int b) {
    return a > b;   // Ascending order, if it is changed to return a>b, it is descending order
}

template<typename T>
// Binary Search
int search_Binary_left(vector<T>v, T value, int low = 0, int high = INT_MAX) {  // search_Binary_left(vector, number to find, left index)

	int vectorSize = v.size()-1;
	high = min(vectorSize, INT_MAX);
	int mid = (low + high) / 2;
	while (low <= high) {
		if (v[mid] == value) {
			return mid;
		} else if (value < v[mid]) {
			high = mid-1;
		} else {
			low = mid + 1;
		}
		mid= (low + high) / 2;
	}

	if (high < 0) {
		high = 0;
	}
	
	return high;
}

template<typename T>
int search_Binary_right(vector<T>v, T value,  int low = 0, int high = INT_MAX) {  // search_Binary_right(vector, number to find, left index)

	int vectorSize = v.size()-1;
	high = min(vectorSize, INT_MAX);
	int mid = (low + high) / 2;
	while (low <= high) {
		if (v[mid] == value) {
			return mid;
		} else if (value < v[mid]) {
			high = mid-1;
		} else {
			low = mid + 1;
		}
		mid= (low + high) / 2;
	}

	if (low > (v.size() - 1)) {
		low = (v.size() - 1);
	}
	
	return low;
}

#endif