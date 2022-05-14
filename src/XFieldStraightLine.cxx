#include <iostream>
#include <iomanip>
#include "XLogger.h"
#include "XFieldStraightLine.h"

using namespace std;

XFieldStraightLine :: XFieldStraightLine()
  : fFlag( -99),
    fNmax(   8),
    fRy0 (-99.),
    fRy1 (-99.),
    fMu_r(1e+3),
    fCurr(1.  ),
    fCart(false)
{
  fPoints = VectorXd :: Zero(4);
}

void XFieldStraightLine :: CartesianCoordinate()
{
  fCart = true;
  Info( "SWITCHED THE FIELD COORDINATE TO CARTESIAN COORDINATE." );
}

void XFieldStraightLine :: SetPoint(const double x, const double y, const double curr)
{
  // push back the points
  fPoints(0) = x;
  fPoints(1) = y;
  fPoints(2) = sqrt(x*x+y*y);
  fPoints(3) = atan2( y, x ); 
  fCurr      = curr;

  // debug
  Debug( setw(15) << setprecision(6) << scientific << x
      << setw(15) << setprecision(6) << scientific << y
      << setw(10) << setprecision(2) << fixed << curr );
}

Vector2d XFieldStraightLine :: GetMagneticField(const double x, const double y)
{
  Vector2d field(0.,0.);

  double r, phi, br, bphi;
  int model = 1000000;

  // select the model
  if      (fRy0<0 && fRy1<0) model = 0;
  else if (fRy0>0 && fRy1<0) model = 1;
  else if (fRy0>0 && fRy1>0) model = 2;

  // coordinate transformation
  r   = sqrt(x*x + y*y);
  phi = atan2(y,x);

  // calculate field
  switch (model) {
    case 0: calc_field_from_wire(r, phi, br, bphi); break;
    case 1: calc_field_infinite_yoke(r, phi, fRy0, fMu_r, br, bphi); break;
    case 2: calc_field_shell_yoke(r, phi, fRy0, fRy1, fMu_r, br, bphi); break;
    default: Fatal( "NO MODEL IS SELECTED FOR FIELD CALCULATION." );
             Fatal( "MODEL ID:" << model );
             break;
  }

  // coordinate transformation for field
  if (fCart) {
    field(0) = cos(phi)*br - sin(phi)*bphi;
    field(1) = sin(phi)*br + cos(phi)*bphi;
  }
  else {
    field(0) = br;
    field(1) = bphi;
  }

  return field;
}

void XFieldStraightLine :: calc_field_from_wire(const double r, const double phi, double& br, double& bphi)
{
  br   = 0.;
  bphi = 0.;

  // for the case r<r'
  if (r<fPoints(2)) {
    for (int j=1; j<fNmax; j++) {
      br  +=  pow(r/fPoints(2), j-1) * sin(j*(fPoints(3)-phi));
      bphi+= -pow(r/fPoints(2), j-1) * cos(j*(fPoints(3)-phi));
    }
  }

  // for the case r>r'
  if (r>fPoints(2)) {
    for (int j=0; j<fNmax; j++) {
      br  += pow(fPoints(2)/r, j+1) * sin(j*(fPoints(3)-phi)); 
      bphi+= pow(fPoints(2)/r, j+1) * cos(j*(fPoints(3)-phi)); 
    }
  }

  br  *= 2e-7*fCurr/fPoints(2);
  bphi*= 2e-7*fCurr/fPoints(2);
}

void XFieldStraightLine :: calc_field_infinite_yoke(const double r, const double phi, const double Rf, const double mu_r, double& br, double& bphi)
{
  br   = 0.; 
  bphi = 0.;

  const double mu_factor = (mu_r-1.) / (mu_r+1.);
  double r_s, phi_s, bbr, bbp;

  r_s   = fPoints(2);
  phi_s = fPoints(3);
  bbr   = 0.;
  bbp   = 0.;

  // for the case r'<r<Rf
  if ( r>r_s && r<Rf ) {
    for (int j=1; j<fNmax; j++) {
      bbr += pow(r_s/r, j+1) * sin(j*(phi_s-phi)) * (1.+mu_factor*pow(r/Rf,2*j));
      bbp += pow(r_s/r, j+1) * cos(j*(phi_s-phi)) * (1.-mu_factor*pow(r/Rf,2*j));
    }

    bbr *= 2e-7*fCurr/r_s;
    bbp *= 2e-7*fCurr/r_s;
    bbp += 2e-7*fCurr/r  ;
  }

  // for the case r<r'
  if ( r < r_s ) {
    for (int j=1; j<fNmax; j++) {
      bbr += pow(r/r_s, j-1) * sin(j*(phi_s-phi)) * (1.+mu_factor*pow(r_s/Rf,2*j));
      bbp += pow(r/r_s, j-1) * cos(j*(phi_s-phi)) * (1.+mu_factor*pow(r_s/Rf,2*j));
    }

    bbr *=  2e-7*fCurr/r_s;
    bbp *= -2e-7*fCurr/r_s;
  }

  // for the case r>Rf
  if ( r > Rf ) {
    for (int j=1; j<fNmax; j++) {
      bbr += pow(r_s/r, j+1) * sin(j*(phi_s-phi));
      bbp += pow(r_s/r, j+1) * cos(j*(phi_s-phi));
    }

    bbr *= 4e-7*mu_r*fCurr / (mu_r+1.) / r_s;
    bbp *= 4e-7*mu_r*fCurr / (mu_r+1.) / r_s;
    bbp += 2e-7*mu_r*fCurr / r;
  }

  br   += bbr;
  bphi += bbp;
}

void XFieldStraightLine :: calc_field_shell_yoke(const double r, const double phi, const double Rf, const double Ra, const double mu_r, double& br, double& bphi)
{
  br   = 0.;
  bphi = 0.;

  const double mu_factor = (mu_r-1.) / (mu_r+1.);
  const double R_factor  = Rf / Ra;
  double term, r_s, phi_s, bbr, bbp;

  r_s   = fPoints(2);
  phi_s = fPoints(3);
  bbr   = 0.;
  bbp   = 0.;

  // for the case r<r'
  if ( r < r_s ) {
    for (int j=1; j<fNmax; j++) {
      term = (1.-pow(R_factor,2*j)) / (1.-pow(mu_factor,2)*pow(R_factor,2*j));
      bbr += pow(r/r_s, j-1) * sin(j*(phi_s-phi)) * (1. + mu_factor*pow(r_s/Rf,2*j)*term);
      bbp += pow(r/r_s, j-1) * cos(j*(phi_s-phi)) * (1. + mu_factor*pow(r_s/Rf,2*j)*term);
    }

    bbr *=  2e-7*fCurr/r_s;
    bbp *= -2e-7*fCurr/r_s;
  }

  // for the case r' < r < Rf
  if ( r>r_s && r<Rf ) {
    for (int j=1; j<fNmax; j++) {
      term = (1.-pow(R_factor,2*j)) / (1.-pow(mu_factor,2)*pow(R_factor,2*j));
      bbr += pow(r_s/r, j+1) * sin(j*(phi_s-phi)) * (1. + mu_factor*pow(r/Rf,2*j)*term);
      bbp += pow(r_s/r, j+1) * cos(j*(phi_s-phi)) * (1. - mu_factor*pow(r/Rf,2*j)*term);
    }

    bbr *= 2e-7 * fCurr / r_s;
    bbp *= 2e-7 * fCurr / r_s;
    bbp += 2e-7 * fCurr / r;
  }

  // for the case Rf < r < Ra
  if ( r>Rf && r<Ra ) {
    for (int j=1; j<fNmax; j++) {
      term  = 1.-pow(mu_factor,2)*pow(R_factor,2*j);
      bbr += pow(r_s/r, j+1) * sin(j*(phi_s-phi)) * (1.-mu_factor*pow(r/Ra,2*j)) / term;
      bbp += pow(r_s/r, j+1) * cos(j*(phi_s-phi)) * (1.+mu_factor*pow(r/Ra,2*j)) / term;
    }

    bbr *= 4e-7*mu_r*fCurr / (r_s*(mu_r+1.));
    bbp *= 4e-7*mu_r*fCurr / (r_s*(mu_r+1.));
    bbp += 2e-7*mu_r*fCurr / r;
  }

  // for the case r > Ra
  if ( r > Ra ) {
    for (int j=1; j<fNmax; j++) {
      term  = 1.-pow(mu_factor,2)*pow(R_factor,2*j);
      bbr += pow(r_s/r, j+1) * sin(j*(phi_s-phi)) / term;
      bbp += pow(r_s/r, j+1) * cos(j*(phi_s-phi)) / term;
    }

    bbr *= 4e-7*mu_r*fCurr / (r_s*pow(mu_r+1,2));
    bbp *= 4e-7*mu_r*fCurr / (r_s*pow(mu_r+1,2));
    bbp += 2e-7*fCurr / r;
  }
   
  br   += bbr;
  bphi += bbp;
}

