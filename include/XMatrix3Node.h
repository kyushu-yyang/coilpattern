#ifndef XMATRIX3NODE_HH
#define XMATRIX3NODE_HH

#include <Eigen/Dense>

using namespace Eigen;

class XMeshHandle2;

/// @file   XMatrix3Node.h
/// @author Ye Yang (QST)
/// @date   04.15.2022

class XMatrix3Node
{
  public:
    /// @brief constructor
    XMatrix3Node();

    /// @brief deconstructor
    ~XMatrix3Node();

    /// @brief setup the order of gauss quad
    void SetOrderOfGauss(const int ngauss);

    /// @brief setup mesh of mandrel
    void SetMeshInfo(XMeshHandle2* handle) { fHandle = handle; }

    /// @brief setup total number of target point
    void SetNumOfTargetMesh(const int num);

    /// @brief append points of target field region
    void SetPointAtTarget(const int i, const double x, const double y, const double z);

    /// @brief return the mesh points at target region
    MatrixXd GetPointAtTarget() { return fPosTar; }

    /// @brief return the response matrix of magnetic field
    void GetFieldMatrix(MatrixXd& bxMat, MatrixXd& byMat, MatrixXd& bzMat);

  protected:
    /// @brief calculate basis vector potential
    void calc_an_vec(const int node_id, const double x, const double y, const double z, Vector3d& an);

    /// @brief calculate basis vector of magnetic field
    void calc_bn_vec(const int node_id, const double x, const double y, const double z, Vector3d& bn);

    /// @brief integration of element area
    Vector3d calc_an_element(const int n_id, XElement3Node* element, Vector3d& rvec);

    Vector3d calc_bn_element(const int n_id, XElement3Node* element, Vector3d& rvec);
    Vector3d calc_bn_element_gauss(const int n_id, XElement3Node* element, Vector3d& rvec);

  private:
    XMeshHandle2* fHandle;
    MatrixXd      fPosTar;
    MatrixXd      fGaussX;
    VectorXd      fGaussW;
};

#endif
