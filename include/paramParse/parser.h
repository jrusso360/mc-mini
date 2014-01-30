#pragma once

#include <vector>
#include <string>

#include "paramParse/paramTree.h"

using namespace std;

class paramParser {
  public:
    paramParser (string);
    ~paramParser ();

    bool push (string);
    bool tryPush (string);
    void pop ();

    void getParamString (string, string&);
    void queryParamString (string, string&, const string);
    void getParamStringVect (string, vector<string>&);
    void queryParamStringVect (string, vector<string>&, const vector<string>);

    void getParamInt (string, int&);
    void queryParamInt (string, int&, const int);
    void getParamIntVect (string, vector<int>&);
    void queryParamIntVect (string, vector<int>&, const vector<int>);

    void getParamDouble (string, double&);
    void queryParamDouble (string, double&, double);
    void getParamDoubleVect (string, vector<double>&);
    void queryParamDoubleVect (string, vector<double>&, const vector<double>);
    
  private:
    paramTree * treeBase;
};

vector<string> inline stringSplit (const string &source, const char *delimiter = " ", bool keepEmpty = false) {
  vector<string> results;

  size_t prev = 0;
  size_t next = 0;

  while ((next = source.find_first_of(delimiter, prev)) != string::npos) {
    if (keepEmpty || (next - prev != 0)) {
      results.push_back(source.substr(prev, next-prev));
    }
    prev = next + 1;
  }
  if (prev < source.size()) {
    results.push_back(source.substr(prev));
  }
  return results;
}
