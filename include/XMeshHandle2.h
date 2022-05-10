#ifndef XMESHHANDLE2_HH
#define XMEHSHANDLE2_HH

#include <Eigen/Dense>
#include <vector>
#include "XVtkFile.h"
#include "XElement3Node.h"

using namespace std;

// @file   XMeshHandle2.h
// @author Y.Yang (QST)
// @date   2022.04.26 (created)

class XMeshHandle2
{
  public:
    /// @brief constructor
    XMeshHandle2();

    /// @brief deconstructor
    ~XMeshHandle2();

    /// @brief setup vtk file handler
    void SetVtkFile(XVtkFile* file);

    /// @brief return the element surrounding the node
    vector<XNode*> GetSurroundingNodes(const int element);

    /// @brief return the node surrounding the element
    vector<XElement3Node*> GetSurroundingElements(const int node);

    /// @brief return the number of nodes
    int GetNumOfNodes() const { return fNode.size(); }

    /// @brief return the node at given index
    XNode* GetNode(const int index) { return fNode.at(index); }

    /// @brief return the number of cells
    int GetNumOfElements() const { return fElement.size(); }

    /// @brief return the cell at given index
    XElement3Node* GetElement(const int index) { return fElement.at(index); }

    /// @brief print out information
    void Print();

  protected:
    /// @brief setup the node vectors
    void setup_nodes(XVtkFile* file);

    /// @brief setup boundary
    void setup_boundaries(XVtkFile* file);

    /// @brief setup element
    void setup_elements(XVtkFile* file);

    /// @breif check the current vector is equal to zero or not
    bool sum_current_vector(const int nid, const vector<XElement3Node*>& elements);

    /// @brief print out the elements surrounding the given node
    void print_elements_sur();

  private:
    vector<XNode*> fNode;
    vector<XElement3Node*> fElement;
};

#endif
