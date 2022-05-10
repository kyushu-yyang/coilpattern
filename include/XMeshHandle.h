#ifndef XMESHHANDLE_HH
#define XMESHHANDLE_HH

#include <Eigen/Dense>
#include <vector>
#include <string>
#include "XElement3Node.h"

using namespace Eigen;
using namespace std;

class XMeshHandle
{
  public:
    /// @brief constructor
    XMeshHandle(const char* filename);

    /// @brief deconstructor
    ~XMeshHandle();

    /// @brief load vtk file
    void Load(const char* filename);

    /// @brief return the number of nodes
    int GetNumOfNodes() const { return fNodes.size(); }

    /// @brief return the number of cells
    int GetNumOfElements() const { return fElements.size(); }

    /// @brief return the node at given index
    XNode* GetNode(const int index) { return fNodes.at(index); }

    /// @brief return the cell at given index
    XElement3Node* GetElement(const int index) { return fElements.at(index); }

    /// @brief search for the elements surrounding the given node
    vector<XElement3Node*> GetSurroundingElements(const int node);

    /// @brief return the node point at the given element
    vector<XNode*> GetSurroundingNodes(const int element);

  protected:
    /// @brief load file line by line
    void load_vtk_file(const char* filename, vector<vector<string>>& table);

  private:
    vector<XNode*> fNodes;
    vector<XElement3Node*> fElements;

};

#endif

