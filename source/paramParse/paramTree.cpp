#include <iostream>
#include <string>
#include <map>

#include "paramParse/parserException.h"
#include "paramParse/paramTree.h"

using namespace std;

paramNode::~paramNode () {
  for (map<string, paramNode *>::iterator it = children.begin(); it != children.end(); ++it) {
    delete (it->second);
  }
}

paramTree::paramTree () {
  rootNode = focusNode = new paramNode;
  focusNode->parent = NULL;
  cerr << "Initialized parameter tree." << endl;
}

paramTree::~paramTree () {
  delete rootNode;
}

void paramTree::moveUp (string key) {
  if (focusNode->children.count (key)) {
    map<string, paramNode *>::iterator iter = focusNode->children.find (key);

    focusNode = iter->second;
  } else {
    throw child_not_found;
  } 
}

void paramTree::moveDown () {
  if (focusNode->parent == NULL) {
    throw -1;
  } else {
    focusNode = focusNode->parent;
  }
}

string paramTree::getParam (string key) {
  if (focusNode->params.count (key)) {
    return focusNode->params[key];
  } else {
    throw param_not_found;
  }
}

void paramTree::addNode (string key) {
  if (focusNode->children.count (key)) {
    throw -1;
  } else {
    paramNode * newNode = new paramNode;
    newNode->parent = focusNode;
    focusNode->children.insert(pair<string, paramNode *> (key, newNode));
  }
}

void paramTree::delNode (string key) {
  if (focusNode->children.count (key)) {
    paramNode * match = focusNode->children.find (key)->second;
    focusNode->children.erase (key);
    delete match;
  } else {
    throw child_not_found;
  }
}

void paramTree::addParam (string key, string param) {
  focusNode->params.insert(pair<string, string> (key, param));
  cerr << "Wrote " << key << " as " << param << endl;
}

void paramTree::delParam (string key) {
  if (focusNode->params.count (key)) {
    focusNode->params.erase (key);
  } else {
    throw param_not_found;
  }
}
