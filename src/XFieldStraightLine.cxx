#include <iostream>
#include <iomanip>
#include "XFieldStraightLine.h"

XFieldStraightLine :: XFieldStraightLine()
  : fFlag( -99),
    fNmax(   8),
    fRy0 (-99.),
    fRy1 (-99.),
    fMu_r(1e+3)
{}

void XFieldStraightLine :: SetPoint(const int idx, const double x, const double y, const double curr)
{
  /// resize the lists
  if (fCurr.size()!=idx+1) {
    fPoints.conservativeResize( idx, 4 );  
    fCurr  .conservativeResize( idx );
  }

  // push back the points
  fPoints(idx, 0) = x;
  fPoints(idx, 1) = y;
  fPoints(idx, 2) = sqrt(x*x+y*y);
  fPoints(idx, 3) = atan2( y, x ); 

  fCurr(idx) = curr;
}

void XFieldStraightLine :: calc_field_from_wire(const double r, const double phi, double& br, double& bphi)
{
  br   = 0.;
  bphi = 0.;

  // sum up the field derived from each wire
  for (int i=0; i<fCurr.size(); i++) {
    // for the case r<r'
    if (r<fPoints(i,2)) {
      for (int j=1; j<fNmax; j++) {
        br  +=  pow(r/fPoints(i,2), j-1) * sin(j*(fPoints(i,3)-phi)) * 2e-7 * fCurr(i) / fPoints(i,2);
        bphi+= -pow(r/fPoints(i,2), j-1) * cos(j*(fPoints(i,3)-phi)) * 2e-7 * fCurr(i) / fPoints(i,2);
      }
    }

    // for the case r>r'
    if (r>fPoints(i,2)) {
      for (int j=0; j<fNmax; j++) {
        br  += pow(fPoints(i,2)/r, j+1) * sin(j*(fPoints(i,3)-phi)) * 2e-7 * fCurr(i) / fPoints(i,2); 
        bphi+= pow(fPoints(i,2)/r, j+1) * cos(j*(fPoints(i,3)-phi)) * 2e-7 * fCurr(i) / fPoints(i,2); 
      }
    }

  }
}

void XFieldStraightLine :: calc_field_infinite_yoke(const double r, const double phi, const double Rf, const double mu_r, double& br, double& bphi)
{
  br   = 0.; 
  bphi = 0.;

  const double mu_factor = (mu_r-1.) / (mu_r+1.);
  double r_s, phi_s, bbr, bbp;

  for (int i=0; i<fCurr.size(); i++) {
    r_s   = fPoints(i,2);
    phi_s = fPoints(i,3);
    bbr   = 0.;
    bbp   = 0.;

    // for the case r'<r<Rf
    if ( r>r_s && r<Rf ) {
      for (int j=1; j<fNmax; j++) {
        bbr += pow(r_s/r, j+1) * sin(j*(phi_s-phi)) * (1.+mu_factor*pow(r/Rf,2*j));
        bbp += pow(r_s/r, j+1) * cos(j*(phi_s-phi)) * (1.-mu_factor*pow(r/Rf,2*j));
      }

      bbr *= 2e-7*fCurr(i)/r_s;
      bbp *= 2e-7*fCurr(i)/r_s;
      bbp += 2e-7*fCurr(i)/r  ;
    }

    // for the case r<r'
    if ( r < r_s ) {
      for (int j=1; j<fNmax; j++) {
        bbr += pow(r/r_s, j-1) * sin(j*(phi_s-phi)) * (1.+mu_factor*pow(r_s/Rf,2*j));
        bbp += pow(r/r_s, j-1) * cos(j*(phi_s-phi)) * (1.+mu_factor*pow(r_s/Rf,2*j));
      }

      bbr *=  2e-7*fCurr(i)/r_s;
      bbp *= -2e-7*fCurr(i)/r_s;
    }

    // for the case r>Rf
    if ( r > Rf ) {
      for (int j=1; j<fNmax; j++) {
        bbr += pow(r_s/r, j+1) * sin(j*(phi_s-phi));
        bbp += pow(r_s/r, j+1) * cos(j*(phi_s-phi));
      }

      bbr *= 4e-7*mu_r*fCurr(i) / (mu_r+1.) / r_s;
      bbp *= 4e-7*mu_r*fCurr(i) / (mu_r+1.) / r_s;
      bbp += 2e-7*mu_r*fCurr(i) / r;
    }

    br   += bbr;
    bphi += bbp;
  }
}

void XFieldStraightLine :: calc_field_shell_yoke(const double r, const double phi, const double Rf, const double Ra, const double mu_r, double& br, double& bphi)
{
  br   = 0.;
  bphi = 0.;

  const double mu_factor = (mu_r-1.) / (mu_r+1.);
  const double R_factor  = Rf / Ra;
  double term, r_s, phi_s, bbr, bbp;

  for (int i=0; i<fCurr.size(); i++) {
    r_s   = fPoints(i,2);
    phi_s = fPoints(i,3);
    bbr   = 0.;
    bbp   = 0.;

    // for the case r<r'
    if ( r < r_s ) {
      for (int j=1; j<fNmax; j++) {
        term = (1.-pow(R_factor,2*j)) / (1.-pow(mu_factor,2)*pow(R_factor,2*j));
        bbr += pow(r/r_s, j-1) * sin(j*(phi_s-phi)) * (1. + mu_factor*pow(r_s/Rf,2*j)*term);
        bbp += pow(r/r_s, j-1) * cos(j*(phi_s-phi)) * (1. + mu_factor*pow(r_s/Rf,2*j)*term);
      }

      bbr *=  2e-7*fCurr(i)/r_s;
      bbp *= -2e-7*fCurr(i)/r_s;
    }

    // for the case r' < r < Rf
    if ( r>r_s && r<Rf ) {
      for (int j=1; j<fNmax; j++) {
        term = (1.-pow(R_factor,2*j)) / (1.-pow(mu_factor,2)*pow(R_factor,2*j));
        bbr += pow(r_s/r, j+1) * sin(j*(phi_s-phi)) * (1. + mu_factor*pow(r/Rf,2*j)*term);
        bbp += pow(r_s/r, j+1) * cos(i*(phi_s-phi)) * (1. - mu_factor*pow(r/Rf,2*j)*term);
      }

      bbr *= 2e-7 * fCurr(i) / r_s;
      bbp *= 2e-7 * fCurr(i) / r_s;
      bbp += 2e-7 * fCurr(i) / r;
    }

    // for the case Rf < r < Ra
    if ( r>Rf && r<Ra ) {
      for (int j=1; j<fNmax; j++) {
        term  = 1.-pow(mu_factor,2)*pow(R_factor,2*j);
        bbr += pow(r_s/r, j+1) * sin(j*(phi_s-phi)) * (1.-mu_factor*pow(r/Ra,2*j)) / term;
        bbp += pow(r_s/r, j+1) * cos(j*(phi_s-phi)) * (1.+mu_factor*pow(r/Ra,2*j)) / term;
      }

      bbr *= 4e-7*mu_r*fCurr(i) / (r_s*(mu_r+1.));
      bbp *= 4e-7*mu_r*fCurr(i) / (r_s*(mu_r+1.));
      bbp += 2e-7*mu_r*fCurr(i) / r;
    }

    // for the case r > Ra
    if ( r > Ra ) {
      for (int j=1; j<fNmax; j++) {
        term  = 1.-pow(mu_factor,2)*pow(R_factor,2*j);
        bbr += pow(r_s/r, j+1) * sin(j*(phi_s-phi)) / term;
        bbp += pow(r_s/r, j+1) * cos(j*(phi_s-phi)) / term;
      }

      bbr *= 4e-7*mu_r*fCurr(i) / (r_s*pow(mu_r+1,2));
      bbp *= 4e-7*mu_r*fCurr(i) / (r_s*pow(mu_r+1,2));
      bbp += 2e-7*fCurr(i) / r;
    }
   
    br   += bbr;
    bphi += bbp;
  }
}

