#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include "XFieldStraightLine.h"
#include "XSolveProblemTSVD.h"
#include "XLogger.h"

using namespace std;
using namespace Eigen; 

void generate_wire(const int layer, const double dLayer, VectorXi& flag, VectorXd& x, VectorXd& y)
{
  // generate an elliptical distribution
  const double a    = 0.075 + dLayer;
  const double b    = 0.055 + dLayer;
  const double eta0 = atanh(b/a);
  const double e    = sqrt(a*a - b*b);
  const double phi0 = 0.;
  const double phi1 = 2*M_PI;
  const int    nphi = 100;
  const double dphi = (phi1-phi0) / (nphi-1);

  const int lenlist = flag.size();

  if (lenlist==0) {
    x    = VectorXd :: Zero(nphi);
    y    = VectorXd :: Zero(nphi);
    flag = VectorXi :: Zero(nphi);
  }
  else {
    x   .conservativeResize( lenlist+nphi ); 
    y   .conservativeResize( lenlist+nphi ); 
    flag.conservativeResize( lenlist+nphi ); 
  }

  double phi, xx, yy;

  Info("APPENDING WIRE AT A: " << setprecision(4) << fixed << a << "[m], "
        << "B: " << setprecision(4) << fixed << b << "[m], "
        << "ETA0:" << setprecision(3) << fixed << eta0 );

  for (int i=lenlist; i<lenlist+nphi; i++) {
    phi = phi0 + dphi*i;
    xx  = e * cosh(eta0) * cos(phi);
    yy  = e * sinh(eta0) * sin(phi);

    x(i) = xx;
    y(i) = yy;
    flag(i) = layer;
  }
}

void generate_voi(VectorXd& x, VectorXd& y, VectorXd& by)
{
  // generate points at Region A
  const double x0_A = -0.019;
  const double x1_A =  0.019;
  const int    nx_A = 10;
  const double dx_A = (x1_A-x0_A) / (nx_A-1);
  const double y0_A = -0.030;
  const double y1_A =  0.030;
  const int    ny_A = 10;
  const double dy_A = (y1_A-y0_A) / (ny_A-1);
  const double bby  = 3.5;

  Info("FIELD AT VOLUME OF INTEREST:" << bby << "T.");
  Info("REGION A: " << setprecision(1) << fixed << (x1_A-x0_A)*1e+3
       << "x" << setprecision(1) << fixed << (y1_A-y0_A)*1e+3 << "mm^2." );

  // generate points at Region B
  const double x0_B = -0.047;
  const double x1_B =  0.047;
  const int    nx_B = 15;
  const double dx_B = (x1_B-x0_B) / (nx_B-1);
  const double y0_B = -0.008;
  const double y1_B =  0.008;
  const int    ny_B = 10;
  const double dy_B = (y1_B-y0_B) / (ny_B-1);

  Info("REGION B: " << setprecision(1) << fixed << (x1_B-x0_B)*1e+3
       << "x" << setprecision(1) << fixed << (y1_B-y0_B)*1e+3 << "mm^2." );

  // generate points
  int npoints = nx_A*ny_A + nx_B*ny_B;
  x = VectorXd :: Zero( npoints );
  y = VectorXd :: Zero( npoints );
  by= VectorXd :: Zero( npoints );
  double xx, yy;
  int cnt = 0;

  for (int i=0; i<nx_A; i++) {
    xx = x0_A + dx_A*i;
    for (int j=0; j<ny_A; j++) {
      yy = y0_A + dy_A*j;
      x (cnt) = xx;
      y (cnt) = yy;
      by(cnt) = bby;
      cnt ++;
    }
  }

  for (int i=0; i<nx_B; i++) {
    xx = x0_B + dx_B*i;
    for (int j=0; j<ny_B; j++) {
      yy = y0_B + dy_B*j;
      x (cnt) = xx;
      y (cnt) = yy;
      by(cnt) = bby;
      cnt ++;
    }
  }
}

void problem()
{
  XFieldStraightLine* wire = new XFieldStraightLine;
  wire->CartesianCoordinate();
  wire->SetMaxOrder(10);

  VectorXd x, y, by, xs, ys;
  VectorXi flag;
  Vector2d bvec;
  
  generate_wire(1, 0.000, flag, xs, ys);
  generate_wire(2, 0.002, flag, xs, ys);
  generate_voi(x, y, by);

  MatrixXd A = MatrixXd :: Zero( x.size(), xs.size() ); 
  cout << A.rows() << "x" << A.cols() << endl;

  for (int i=0; i<x.size(); i++) {
    for (int j=0; j<xs.size(); j++) {
      wire->SetPoint( xs(j), ys(j) );
      bvec = wire->GetMagneticField( x(i), y(i) );
      A(i,j) = bvec(1); 
    }
  }

  // solve the problem using TSVD
  XSolveProblemTSVD* sol = new XSolveProblemTSVD;
  sol->SetVolumeOfInterest(x, y, by);
  sol->SetMaxNumber(20);
  sol->Solve(A);
  sol->SetSourcePosition(flag, xs, ys);
  sol->Write("res_ellip_dipole.txt");
}

int main(int argc, char** argv)
{
  // setup logger
  XLogger* logger = XLogger::GetInstance(); 
  logger->Start(XLogger::DEBUG, "one_dimen_debug.log");

  // solve the one-dimensional problem
  try {
    problem();
  }
  catch (invalid_argument& except) {
    cerr << except.what() << endl;
    return EXIT_FAILURE;
  }

  logger->Stop();

  return EXIT_SUCCESS;
}
