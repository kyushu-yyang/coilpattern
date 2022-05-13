#ifndef XELEMENT_HH
#define XELEMENT_HH

#include <Eigen/Dense>
#include <vector>
#include "XNode.h"

using namespace Eigen;
using namespace std;

/// @file   XElement.h
/// @author Ye Yang (QST)
/// @date   05.13.2022

/// base class
class XElement
{
  public:
    /// @brief constructor
    XElement();

    /// @brief deconstructor
    virtual ~XElement();

    /// @brief setup element id number
    void SetId(const int id) { fId = id; }

    /// @brief return the element id number 
    int  GetId() const { return fId; }

    /// @brief setup surface id number
    void SetSurfaceId(const int id);

    /// @brief return the id of surface
    int  GetSurfaceId() const { return fSurfId; }

    /// @brief setup point of node
    void SetNode(const int index, XNode* node);

    /// @brief return the node
    vector<XNode*> GetNodes() { return fNode; }

    /// @brief return the node at the given node id
    void GetNode(const int id, XNode* node);

    /// @brief return the node at the given index
    XNode* GetNode(const int index);

    /// @brief check whether the element is valid or not
    bool Is_Valid() const { return fValid; }

  protected:
    vector<XNode*> fNode;

  private:
    int  fId;
    int  fSurfId;
    bool fValid;
};

#endif
