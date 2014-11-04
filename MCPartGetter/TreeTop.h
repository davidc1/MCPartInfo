/**
 * \file TreeTop.h
 *
 * \ingroup MCInfo
 * 
 * @author David Caratelli
 */

/** \addtogroup MCPartGetter

    @{*/
#ifndef TREETOP_H
#define TREETOP_H

#include <iostream>
#include "TreeNode.h"

/**
   \class TreeTop
   User defined class TreeTop ... these comments are used to generate
   doxygen documentation!
 */
class TreeTop : public TreeNode {
  
public:

  TreeTop() : TreeNode() {};

  TreeTop(int index) : TreeNode() { _index = index; };

  
private:
  
};

#endif
/** @} */ // end of doxygen group 

