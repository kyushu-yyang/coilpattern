#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <Eigen/Dense>
#include "XFieldStraightLine.h"
#include "XLogger.h"

using namespace std;
using namespace Eigen; 

const double RREF = 0.030;

void load_file(const string filename, vector<double>& x, vector<double>& y, vector<double>& I)
{
  ifstream inputfile( filename );

  // check file exists or not
  if ( !inputfile.is_open() ) {
    cerr << "ERROR: Could not open the file - " << filename << endl;
    Fatal("CANNOT OPEN THE DATA FILE: " << filename << ".");
    throw invalid_argument("could not open the file.");
  }

  if (x.size()!=0) {
    x.clear();
    y.clear();
    I.clear();
  }

  string eachline;
  int cnt = 0;

  while (getline(inputfile, eachline)) {
    istringstream iss(eachline);
    int    flag, curr;
    double xx, yy;

    if (!(iss >> xx >> yy >> curr))
      break;

    x.push_back( xx );
    y.push_back( yy );
    I.push_back( curr);

    cnt ++;
  }

  cout << "TOTAL NUMBER OF LINE CURRENT:" << cnt << endl;
  cout << "FINISHED" << endl;
}

void generate_voi(VectorXd& x, VectorXd& y, VectorXd& bx, VectorXd& by)
{
  const double theta0 = 0.;
  const double theta1 = 2.*M_PI;
  const int    ntheta = 51;
  const double dtheta = (theta1-theta0)/(ntheta-1);
  double       theta  = 0.;

  x = VectorXd :: Zero(ntheta);
  y = VectorXd :: Zero(ntheta);
  bx= VectorXd :: Zero(ntheta);
  by= VectorXd :: Zero(ntheta);

  for (int i=0; i<ntheta; i++) {
    theta = theta0 + i*dtheta;
    x(i)  = RREF * cos(theta);
    y(i)  = RREF * sin(theta);
  }
}

void generate_point_at_source(const vector<double>& xs, const vector<double>& ys, const vector<double>& Is,
                              vector<double>& x_s, vector<double>& y_s, vector<double>& bx_s, vector<double>& by_s)
{
  x_s .clear();
  y_s .clear();
  bx_s.clear();
  by_s.clear();

  for (int i=0; i<xs.size(); i++) {
    if (Is.at(i)!=0.) {
      x_s.push_back( xs.at(i) );
      y_s.push_back( ys.at(i) );
      bx_s.push_back( 0. );
      by_s.push_back( 0. );
    }
  }
}

void get_multipole(const VectorXd& bx, const VectorXd& by, int n, double& an, double& bn)
{
  const double theta0 = 0.;
  const double theta1 = 2.*M_PI;
  const int    ntheta = 51;
  const double dtheta = (theta1-theta0)/(ntheta-1);
  double       t0, t1;
  complex<double> b0, b1, exp0, exp1;
  complex<double> cn(0.,0.);

  for (int i=1; i<ntheta; i++) {
    t0 = theta0 + (i-1)*dtheta;
    t1 = theta0 + i*dtheta;

    b0.real(by(i-1)); b0.imag(bx(i-1)); 
    b1.real(by(i  )); b1.imag(bx(i  ));

    exp0.real( cos((n-1)*t0) );  exp0.imag( -sin((n-1)*t0) );
    exp1.real( cos((n-1)*t1) );  exp1.imag( -sin((n-1)*t1) );

    cn += (b0*exp0 + b1*exp1) * dtheta * 0.5;
  }

  cn /= 2*M_PI;

  an = cn.imag();
  bn = cn.real();
}

double get_heat_gen(const double bx, const double by)
{
  const double B0     = 0.34;
  const double Jc0    = 2.0e+10;
  const double B1     = 0.0;
  const double lambda = 1.0/3.83;
  const double rfila  = 1.2e-6;
  const double area   = M_PI*pow(1.0e-3,2)/4;
  const double B2     = sqrt(bx*bx + by*by);

  double Qhys = 8*rfila*Jc0*B0*log((B2+B0)/(B1+B0))*lambda*area/(3*M_PI);

  return Qhys;
}

void problem(const string filename)
{
  XFieldStraightLine* wire = new XFieldStraightLine;
  wire->CartesianCoordinate();
  wire->SetMaxOrder(5);
  wire->SetRelativePermeability(1500);
  wire->SetIronYokeRadius(0.125, -99.0);

  Vector3d bvec;
  VectorXd x, y, bx, by;
  vector<double> xs, ys, Is, x_s, y_s, by_s, bx_s;

  load_file( filename, xs, ys, Is ); 
  generate_voi( x, y, bx, by );
  generate_point_at_source( xs, ys, Is, x_s, y_s, bx_s, by_s );

  // 1/4 model
  for (int i=0; i<x.size(); i++) {
    for (int j=0; j<xs.size(); j++) {
      if (Is.at(j)!=0.) {
        wire->SetPoint( xs.at(j),  ys.at(j),-Is.at(j) );
        bvec = wire->GetMagneticField( x(i), y(i) );
        bx(i) += bvec(0);
        by(i) += bvec(1);

        // flip ymin
        wire->SetPoint( xs.at(j), -ys.at(j),-Is.at(j) );
        bvec = wire->GetMagneticField( x(i), y(i) );
        bx(i) += bvec(0);
        by(i) += bvec(1);

        // flip xmin
        wire->SetPoint(-xs.at(j), ys.at(j), Is.at(j) );
        bvec = wire->GetMagneticField( x(i), y(i) );
        bx(i) += bvec(0);
        by(i) += bvec(1);

        // flip ymin, xmin
        wire->SetPoint(-xs.at(j), -ys.at(j), Is.at(j) );
        bvec = wire->GetMagneticField( x(i), y(i) );
        bx(i) += bvec(0);
        by(i) += bvec(1);
      }
    }
  }

  ofstream outfile("results.txt");
  outfile << "MULTIPOLES:" << endl;

  // calculate multipoles
  const int nmax = 11;
  double an, bn;

  outfile << setw(5 ) << fixed << "N"
          << setw(15) << fixed << "AN[T]" << setw(15) << fixed << "BN[T]" << "\n";
  for (int i=1; i<11; i++) {
    get_multipole( bx, by, i, an, bn );
    outfile << setw(5) << fixed << i << setw(15) << setprecision(7) << scientific << an
            << setw(15) << setprecision(7) << scientific << bn << "\n";
  }

  // calculate field, force, energy and inductance
  double energy = 0.;
  double induct = 0.;
  double Inom;

  for (int i=0; i<x_s.size(); i++) {
    for (int j=0; j<xs.size(); j++) {
      if (Is.at(j)!=0.) {
        wire->SetPoint( xs.at(j),  ys.at(j), Is.at(j) );
        bvec = wire->GetMagneticField( x_s.at(i), y_s.at(i) );
        bx_s.at(i) += bvec(0);
        by_s.at(i) += bvec(1);
        induct += abs(bvec(2));

        // flip ymin
        wire->SetPoint( xs.at(j), -ys.at(j), Is.at(j) );
        bvec = wire->GetMagneticField( x_s.at(i), y_s.at(i) );
        bx_s.at(i) += bvec(0);
        by_s.at(i) += bvec(1);
        induct  += abs(bvec(2));

        // flip xmin
        wire->SetPoint(-xs.at(j), ys.at(j),-Is.at(j) );
        bvec = wire->GetMagneticField( x_s.at(i), y_s.at(i) );
        bx_s.at(i) += bvec(0);
        by_s.at(i) += bvec(1);
        induct  += abs(bvec(2));

        // flip ymin, xmin
        wire->SetPoint(-xs.at(j), -ys.at(j),-Is.at(j) );
        bvec = wire->GetMagneticField( x_s.at(i), y_s.at(i) );
        bx_s.at(i) += bvec(0);
        by_s.at(i) += bvec(1);
        induct  += abs(bvec(2));

        Inom = Is.at(i);
      }
    }
  }

  induct /= Inom;
  energy = 0.5*induct*Inom*Inom;
  outfile << "L[H/m]:" << setw(15) << setprecision(6) << scientific << induct << "\n"; 
  outfile << "E[J/m]:" << setw(15) << setprecision(6) << scientific << energy << "\n"; 

  outfile << setw(15) << "X[m]" << setw(15) << "Y[m]" << setw(15) << "BX[T]" << setw(15) << "BY[T]"
          << setw(15) << "FX[N/m]" << setw(15) << "FY[N/m]" << setw(15) << "Qhys[J/m]" << "\n"; 

  for (int i=0; i<x_s.size(); i++) 
    outfile << setw(15) << setprecision(6) << scientific << x_s.at(i) 
            << setw(15) << setprecision(6) << scientific << y_s.at(i) 
            << setw(15) << setprecision(6) << scientific << bx_s.at(i)
            << setw(15) << setprecision(6) << scientific << by_s.at(i)
            << setw(15) << setprecision(6) << scientific << -Inom*by_s.at(i)
            << setw(15) << setprecision(6) << scientific <<  Inom*bx_s.at(i)
            << setw(15) << setprecision(6) << scientific << get_heat_gen(bx_s.at(i), by_s.at(i)) << endl;
  outfile.close();
}

int main(int argc, char** argv)
{
  const string filename = argv[1]; 

  // solve the one-dimensional problem
  try {
    problem(filename);
  }
  catch (invalid_argument& except) {
    cerr << except.what() << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
