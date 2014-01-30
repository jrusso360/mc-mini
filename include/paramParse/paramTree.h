#pragma once

#include <string>
#include <map>

using namespace std;

class paramNode {
  public:
    ~paramNode();
    paramNode * parent;
    map <string, paramNode *> children;
    map <string, string> params;
};

class paramTree {
  public:
    paramTree ();
    ~paramTree ();

    void moveUp (string);
    void moveDown ();
    string getParam (string);
    void addNode (string);
    void delNode (string);
    void addParam (string, string);
    void delParam (string);
  private:
    paramNode * rootNode;
    paramNode * focusNode;
};
