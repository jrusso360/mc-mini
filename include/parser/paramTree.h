#pragma once

#include <string>
#include <map>

using namespace std;

class ParamNode {
  public:
    ~ParamNode();
    ParamNode * parent;
    map <string, ParamNode *> children;
    map <string, string> params;
};

class ParamTree {
  public:
    ParamTree ();
    ~ParamTree ();

    void moveUp (string);
    void moveDown ();
    string getParam (string);
    void addNode (string);
    void delNode (string);
    void addParam (string, string);
    void delParam (string);
  private:
    ParamNode * rootNode;
    ParamNode * focusNode;
};
