#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/math/special_functions/legendre.hpp>
#include "XTargetSphere.h"

XTargetSphere :: XTargetSphere()
  : fPhi0  (0.), fPhi1  (2*M_PI),
    fTheta0(0.), fTheta1(M_PI  ),
    fR0    (0.), fR1    (1.0   ),
    fNphi  (10), fNtheta(10    ),
    fNr    (1 ), fRref  (0.025 )
{}

void XTargetSphere :: SetPhi(const double phi0, const double phi1, const int nphi)
{
  fPhi0 = phi0; 
  fPhi1 = phi1;
  fNphi = nphi;

  if (fPhi0<0     ) fPhi0 = 0.;
  if (fPhi1>2*M_PI) fPhi1 = 2.*M_PI;
}

void XTargetSphere :: SetTheta(const double theta0, const double theta1, const int ntheta)
{
  fTheta0 = theta0;
  fTheta1 = theta1;
  fNtheta = ntheta;

  if (fTheta0<0   ) fTheta0 = 0.;
  if (fTheta1>M_PI) fTheta1 = M_PI;
}

void XTargetSphere :: SetRadius(const double r0, const double r1, const int nr)
{
  fR0 = r0;
  fR1 = r1;
  fNr = nr;

  if (fR0<0) fR0 = 0.;
}

bool XTargetSphere :: Is_Harmonics_Exist(const int n, const int m)
{
  bool is_exist = false;

  for (int i=0; i<fN.size(); i++) {
    if (fN.at(i)==n && fM.at(i)==m) { is_exist = true; break; }
  }

  return is_exist;
}

void XTargetSphere :: SetHarmonics(const double rref, const int n, const int m, const double Anm, const double Bnm)
{
  fRref = rref;

  // check harmonics exists or not
  if (!Is_Harmonics_Exist(n,m)) {
    fN.push_back( n );
    fM.push_back( m );
    fAnm.push_back( Anm );
    fBnm.push_back( Bnm );
  }
}

MatrixXd XTargetSphere :: GetXYZ()
{
  MatrixXd mat = MatrixXd :: Zero(fNr*fNphi*fNtheta, 3);

  double x, y, z, r, phi, theta;
  int cnt{0};

  const double dphi   = (fPhi1-fPhi0) / (fNphi-1);
  const double dtheta = (fTheta1-fTheta0) / (fNtheta-1);
  const double dr     = (fR1-fR0) / (fNr-1);

  for (int i=0; i<fNr; i++) {
    r = fR0 + i*dr; 
    for (int j=0; j<fNtheta; j++) {
      theta = fTheta0 + j*dtheta;
      for (int k=0; k<fNphi; k++) {
        phi = fPhi0 + k*dphi;
        x   = r * sin(theta) * cos(phi);
        y   = r * sin(theta) * sin(phi);
        z   = r * cos(theta);

        mat(cnt, 0) = x;
        mat(cnt, 1) = y;
        mat(cnt, 2) = z;

        cnt ++;
      }
    }
  }

  return mat;
}

VectorXd XTargetSphere :: GetField()
{
  VectorXd vec = VectorXd :: Zero(fNr*fNphi*fNtheta);

  double r, phi, theta, bz;
  int cnt{0};

  const double dphi   = (fPhi1-fPhi0) / (fNphi-1);
  const double dtheta = (fTheta1-fTheta0) / (fNtheta-1);
  const double dr     = (fR1-fR0) / (fNr-1);

  for (int i=0; i<fNr; i++) {
    r = fR0 + i*dr; 
    for (int j=0; j<fNtheta; j++) {
      theta = fTheta0 + j*dtheta;
      for (int k=0; k<fNphi; k++) {
        phi = fPhi0 + k*dphi;
        bz  = get_magnetic_field(r, theta, phi);
        vec(cnt) = bz;
        cnt ++;
      }
    }
  }

  return vec;
}

double XTargetSphere :: get_magnetic_field(const double r, const double theta, const double phi)
{
  const double u = cos(theta);

  double pnm, anm, bnm;
  int    n, m;
  double bz = 0.;

  for (int i=0; i<fN.size(); i++) {
    n   = fN.at(i);
    m   = fM.at(i);
    anm = fAnm.at(i);
    bnm = fBnm.at(i);
    pnm = boost::math::legendre_p(n, m, u);    

    bz += pow(r/fRref,n) * pnm * (anm*cos(m*phi)+bnm*sin(m*phi));
  }

  return bz;
}
