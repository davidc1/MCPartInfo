#ifndef TREE_CXX
#define TREE_CXX

#include "TreeNode.h"

void TreeNode::countNodes(int& count){

  count += 1;
  for (size_t i=0; i < _children.size(); i++)
    _children.at(i).countNodes(count);

  return;
}

#endif
