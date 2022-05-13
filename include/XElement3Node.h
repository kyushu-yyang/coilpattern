#ifndef XELEMENT3NODE_HH
#define XELEMENT3NODE_HH

#include <Eigen/Dense>
#include <vector>
#include "XNode.h" 

using namespace Eigen;
using namespace std;

class XElement3Node
{
  public:
    /// @brief constructor
    XElement3Node();

    /// @brief deconstructor
    ~XElement3Node();

    /// @brief return the number of node
    int GetNumOfNodes() { return 3; }

    /// @brief setup the id of element
    void SetId(const int id) { fId = id; }

    /// @brief return the id of element
    int  GetId() const { return fId; }

    /// @brief setup surface id
    void SetSurfaceId(const int id) { fSurfId = id; }

    /// @brief return surface id
    int  GetSurfaceId() const { return fSurfId; }

    /// @brief setup the node point
    void SetNode(const int id, XNode* pts) { fNode.at(id) = pts; }

    /// @brief return the node point
    XNode* GetNode(const int index) { return fNode.at(index); }

    /// @brief return the node index opposite to the line of boundary
    int GetBoundaryNodeIndex();

    /// @brief return the node point at given node id
    void GetNode(const int id, XNode* node);

    /// @brief return the node id in this element
    Vector3i GetNodeId();

    /// @brief calculate the shape function at Cartesian coordinates
    Vector3d GetShapeFunction(const double xi, const double eta);

    /// @brief calculate the position vector at given local point
    Vector3d GetPositionVector(const double xi, const double eta);

    /// @brief calculate the area of triangle
    double GetArea();

    /// @brief calculate current basis vector
    MatrixXd GetCurrentBasisVector();
    MatrixXd GetLoopBasisVector();

    /// @breif return current basis vector at given node
    Vector3d GetCurrentBasisVector(const int node);
    Vector3d GetLoopBasisVector(const int node);

  private:
    int fId;
    int fSurfId;
    vector<XNode*> fNode;
};

#endif
