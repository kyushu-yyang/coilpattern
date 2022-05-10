#ifndef XELEMENT3NODE_HH
#define XELEMENT3NODE_HH

#include <Eigen/Dense>
#include <vector>

using namespace Eigen;
using namespace std;

class XNode
{
  public:
    XNode() { fXYZ(0)=0.; fXYZ(1)=0.; fXYZ(2)=0.; }
    ~XNode() {}
    void SetId(const int num) { id = num; }
    int  GetId() const { return id; }
    void   SetPoint(const double xx, const double yy, const double zz) { fXYZ(0)=xx; fXYZ(1)=yy; fXYZ(2)=zz; }
    void   SetX(const double xx) { fXYZ(0)=xx; }
    void   SetY(const double yy) { fXYZ(1)=yy; }
    void   SetZ(const double zz) { fXYZ(2)=zz; }
    double GetX() const { return fXYZ(0); }
    double GetY() const { return fXYZ(1); } 
    double GetZ() const { return fXYZ(2); }
    Vector3d GetPoint() { return fXYZ; }
    Vector3d GetPositionVector(const double x, const double y, const double z);
    double   GetDistance(const double x, const double y, const double z);
    void   AtBoundary(const bool bc) { fNodeAtBC = bc; }
    bool   Is_AtBoundary() { return fNodeAtBC; }

  private:
    int    id{-99};
    Vector3d fXYZ;
    bool   fNodeAtBC{false};
};

/**************************************************/
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
    double fId;
    vector<XNode*> fNode;
};

#endif
