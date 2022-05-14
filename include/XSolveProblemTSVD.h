#ifndef XSOLVEPROBLEMTSVD_HH
#define XSOLVEPROBLEMTSVD_HH

#include <string>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

/// @file   XSolveProblemTSVD.h
/// @author Ye Yang (QST)
/// @date   05.14.2022

class XSolveProblemTSVD
{
  public:
    /// @brief constructor
    XSolveProblemTSVD();

    /// @brief deconstructor
    virtual ~XSolveProblemTSVD() {}

    /// @brief setup the maximum truncated number
    void SetMaxNumber(const int nmax) { fNmax = nmax; }

    /// @brief setup target field region
    void SetVolumeOfInterest(const VectorXd& x, const VectorXd& y, const VectorXd& b);

    /// @breif solve the problem using TSVD
    void Solve(const MatrixXd& A);

    /// @brief setup the position and flags
    void SetSourcePosition(const VectorXi flag, const VectorXd& xs, const VectorXd& ys);

    /// @brief write out the data to text file
    void Write(const string filename);


  protected:
    /// @brief singular value decomposition
    void truncated_SVD(const MatrixXd& A, const VectorXd& b, VectorXd& x);


  private:
    int      fNmax;
    MatrixXd fPos;
    VectorXd fBtar;
    VectorXd fPtg;
    MatrixXd fIvec;
    MatrixXd fBrec;
    VectorXi fFlag;
    MatrixXd fPosSrc;
};

#endif
