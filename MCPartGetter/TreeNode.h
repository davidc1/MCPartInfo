/**
 * \file Tree.h
 *
 * \ingroup MCPartGetter
 * 
 * \brief Class def header for a class Tree
 *
 * @author David Caratelli
 */

/** \addtogroup MCInfo

    @{*/
#ifndef TREE_H
#define TREE_H

#include <iostream>
#include <vector>

/**
   \class Tree
   User defined class Tree ... these comments are used to generate
   doxygen documentation!
 */
class TreeNode{

public:

  /// Default Constructor
  TreeNode(){ _children.clear(); };

  /// constructor to use
  TreeNode(int index){ _index = index; _children.clear(); };

  /// Default destructor
  virtual ~TreeNode(){};

  /// Set if Primary
  void SetPrimary(bool on) { _isPrimary = on; }

  /// is primary?
  bool isPrimary() { return _isPrimary; }

  /// Add a child to this noe
  void AddChild(TreeNode tn) { _children.push_back(tn); };

  /// Get children of this node
  std::vector<TreeNode> getChildren() { return _children; }

  /// Count number of nodes
  void countNodes(int& count);

  /// Get node index
  int getNodeIndex() { return _index; }

  /// Get closest brother (parent's child that is next in the list)
  //TreeNode* getClosestBrother();

  /// Get closest uncle (paren't sibling that is next in the list)
  //TreeNode* getClosestUncle();

  /// Is this node a leaf (i.e. the end of the line...no sub-nodes?)
  bool isLeaf() { if (_children.size() == 0) { return true; } return false; }

  void SetParentId(int id) { _parentID = id; }

  void SetAncestorId(int id) { _ancestorID = id; }

  int getParentId() { return _parentID; }

  int getAncestorId() { return _ancestorID; }


protected:

  //when using this class for MCparticles index should be TrackId
  int _index;

  bool _isPrimary;

  int _parentID;
  
  int _ancestorID;

  std::vector<TreeNode> _children;

};

#endif
/** @} */ // end of doxygen group 

