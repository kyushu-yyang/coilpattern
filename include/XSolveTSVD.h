#ifndef XSOLVETSVD_HH
#define XSOLVETSVD_HH

#include <string>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class XMeshHandle2;

class XSolveTSVD
{
  public:
    /// @brief constructor
    XSolveTSVD();

    /// @brief deconstructor
    ~XSolveTSVD();

    /// @brief setup mesh handler
    void SetMeshHandle(XMeshHandle2* hand) { fHandle = hand; }

    /// @brief setup the truncated mode
    void SetTruncatedNumber(const int nt) { fNumTruncated = nt; }

    /// @brief save result or not
    void SetSaveFile(const string filename);

    /// @brief solve the problem of Ax = b
    void Solve(const MatrixXd& A, const VectorXd& b, VectorXd& x);

  protected:
    /// @brief sigular value decomposition
    void truncated_svd(const MatrixXd& A, const VectorXd& b, VectorXd& x);

    /// @brief save stream function distribution to file
    void save_each_mode(const int mode, const VectorXd& vec); 

    /// @brief reconstruction of magntic field
    double reconstruction(const MatrixXd& A, const VectorXd& b, const VectorXd& x, const bool print);
    
  private:
    XMeshHandle2* fHandle;
    int fNumTruncated;
    string fFilename;
    bool fSave;
};

#endif
