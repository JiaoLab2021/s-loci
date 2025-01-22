#ifndef strip_split_join_hpp
#define strip_split_join_hpp
#include <string>
#include <vector>

using namespace std;


/*
	Remove special characters from both ends of a string
	str -> string
	ch -> special characters
*/
string strip(const string & str, char ch=' ') {
	int i = 0;
	while (str[i] == ch)  // The number of characters in the head is i
		i++;
	int j = str.size() - 1;
	while (str[j] == ch)  // The number of characters at the end is str.size() - 1 - j
		j--;		
	return str.substr(i, j+1-i);
}


/*
	Split the string
	str -> string
	delim -> split character
*/
vector<string> split(const string & str, const string & delim) {
	vector<string> res;  // Store the split substrings in a vector
	if("" == str) return  res;  // If empty, returns an empty vector
	
	string strs = str + delim;  // Expand the string to facilitate retrieval of the last delimited string
	size_t pos;
	size_t size = strs.size();
 
	for (int i = 0; i < size; ++i) {
		pos = strs.find(delim, i); // pos is the position where the separator first appears, and the string from i to before pos is the separated string
		if( pos < size) {  // If found, if no separator is found, pos is string::npos
			string s = strs.substr(i, pos - i);  // Substring starting from i and length pos-i
			res.push_back(s);  // The string cut between two consecutive spaces is an empty string. There is no judgment here whether s is empty, so there is an output of empty characters in the final result.
			i = pos + delim.size() - 1;
		}
	}
	return res;	
}


/*
	String join
	val -> container to be joined
	delim -> characters to be added between container elements
*/
string join(vector<int> & val, string delim) {
    std::string str;
	int vecSize = val.size();
	int index = 0;
	for (auto iter : val) {
		str += to_string(iter);
		
		if (index != vecSize-1)
        {
            str += delim;
        }
		index++;
	}
    return str;
}


/*
	String join
	val -> container to be joined
	delim -> characters to be added between container elements
*/
string join(vector<long int> & val, string delim) {
    std::string str;
	int vecSize = val.size();
	int index = 0;
	for (auto iter : val) {
		str += to_string(iter);
		
		if (index != vecSize-1) {
            str += delim;
        }
		index++;
	}
    return str;
}


/*
	String join
	val -> container to be joined
	delim -> characters to be added between container elements
*/
string join(vector<float> & val, string delim) {
    std::string str;
	int vecSize = val.size();
	int index = 0;
	for (auto iter : val) {
		str += to_string(iter);
		
		if (index != vecSize-1) {
            str += delim;
        }
		index++;
	}
    return str;
}


/*
	String join
	val -> container to be joined
	delim -> characters to be added between container elements
*/
string join(vector<string> & val, string delim) {
    std::string str;
	int vecSize = val.size();
	int index = 0;
	for (auto iter : val) {
		str += iter;
		
		if (index != vecSize-1) {
            str += delim;
        }
		index++;
	}
    return str;
}

#endif