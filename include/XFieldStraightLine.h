#ifndef XFIELDSTRAIGHTLINE_HH
#define XFIELDSTRAIGHTLINE_HH

#include <Eigen/Dense>

using namespace Eigen;

/// @file   XFieldStaightLine.h
/// @author Y.Yang (QST)
/// @date   2022.05.11

class XFieldStraightLine
{
  public:
    /// @brief constructor
    XFieldStraightLine();

    /// @brief deconstructor
    ~XFieldStraightLine() {}

    /// @brief setup flag
    void SetFlag(const int flag) { fFlag = flag; }

    /// @brief return the flag
    int  GetFlag() const { return fFlag; }

    /// @brief setup the limitation of the order for the sum up of field
    void SetMaxOrder(const int nmax) { fNmax = nmax; }

    /// @brief setup the radius of iron yoke
    void SetIronYokeRadius(const double r0, const double r1);

    /// @brief setup the relative permeability for iron yoke
    void SetRelativePermeability(const double mu_r) { fMu_r = mu_r; }

    /// @brief setup the point in two dimension
    void SetPoint(const int idx, const double x, const double y, const double curr=1.);

  protected:
    /// @brief calculate the two-dimensional magnetic field derived from strainght line conductor
    void calc_field_from_wire(const double r, const double phi, double& br, double& bphi);

    /// @brief calculate the two-dimensional magnetic field using infinite iron yoke model
    void calc_field_infinite_yoke(const double r, const double phi, const double Rf, const double mu_r, double& br, double& bphi);

    /// @brief calculate the two-dimensional magnetic field using shell iron yoke model
    void calc_field_shell_yoke(const double r, const double phi, const double Rf, const double Ra, const double mu_r, double& br, double& bphi);

  private:
    int      fFlag;
    int      fNmax;
    double   fRy0;
    double   fRy1;
    double   fMu_r;
    MatrixXd fPoints;
    VectorXd fCurr;

};

#endif
