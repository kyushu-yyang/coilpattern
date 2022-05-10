#ifndef XTARGETSPHERE_HH
#define XTARGETSPHERE_HH

#include <Eigen/Dense>
#include <vector>

using namespace Eigen;
using namespace std;

class XTargetSphere
{
  public:
    /// @brief constructor
    XTargetSphere();

    /// @brief deconstructor
    ~XTargetSphere() {}

    /// @brief setup mesh and range along phi direction
    void SetPhi(const double phi0, const double phi1, const int nphi);

    /// @brief setup mesh and range along theta direction
    void SetTheta(const double theta0, const double theta1, const int ntheta);

    /// @brief setup mesh and range along radial direction
    void SetRadius(const double r0, const double r1, const int nr);
    
    /// @brief setup harmonics of magnetic field
    void SetHarmonics(const double rref, const int n, const int m, const double Anm, const double Bnm);

    /// @brief check the harmonics exists or not
    bool Is_Harmonics_Exist(const int n, const int m);

    /// @brief return the points on sphere
    MatrixXd GetXYZ(); 

    /// @brief return the vector containing the target field
    VectorXd GetField();
 
  protected:
    /// @brief calculate magnetic field from harmonics
    double get_magnetic_field(const double r, const double theta, const double phi);

  private:
    double fPhi0;
    double fPhi1;
    double fTheta0;
    double fTheta1;
    double fR0;
    double fR1;
    int    fNphi;
    int    fNtheta;
    int    fNr;
    double fRref;
    vector<int> fN;
    vector<int> fM;
    vector<double> fAnm;
    vector<double> fBnm;
};

#endif
