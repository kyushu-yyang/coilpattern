#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include "XLogger.h"
#include "XElement.h"

XElement :: XElement()
  : fId    (-9999),
    fSurfId(-9999),
    fValid (false)
{
  fNode.push_back( nullptr );
  fNode.push_back( nullptr );
  fNode.push_back( nullptr );
}

XElement :: ~XElement()
{
  if (fNode.size()!=0) {
    for (auto eachnode : fNode)
      delete eachnode;
    fNode.clear();
  }
}

void XElement :: SetSurfaceId(const int id)
{
  fSurfId = id;
  fValid  = true;
}

void XElement :: SetNode(const int index, XNode* node)
{
  if (index>fNode.size()-1)
    fNode.resize(index+1);

  fNode.at(index) = node;
}

void XElement :: GetNode(const int id, XNode* node)
{
  bool is_node = false;

  for (int i=0; i<fNode.size(); i++) {
    if (id == fNode.at(i)->GetId()) {
      is_node = true;
      node = fNode.at(i);
      break;
    }
  }

  if (!is_node) {
    Fatal("THE INPUTED NODE ID: " << id << " IS OUT OF RANGE.");
    throw out_of_range( "The inputed node id does not exist." );
  }
}

XNode* XElement :: GetNode(const int index)
{
  // check the input
  if (index<0 || index>=fNode.size()) {
    Fatal("THE INPUTED NODE INDEX: " << index << " IS OUT OF RANGE.");
    throw out_of_range("The inputed node index is out of range. please check the index");
  }

  return fNode.at(index);
}
