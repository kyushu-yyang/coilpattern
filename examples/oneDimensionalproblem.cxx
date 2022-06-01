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
  const double a    = 0.080 + dLayer;
  const double b    = 0.060 + dLayer;
  const double eta0 = atanh(b/a);
  const double e    = sqrt(a*a - b*b);
  const double phi0 = -M_PI;
  const double phi1 =  M_PI;
  const int    nphi = 201;
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
  const int    nx_A = 15;
  const double dx_A = (x1_A-x0_A) / (nx_A-1);
  const double y0_A = -0.030;
  const double y1_A =  0.030;
  const int    ny_A = 19;
  const double dy_A = (y1_A-y0_A) / (ny_A-1);
  const double bby  = 3.5;

  Info("FIELD AT VOLUME OF INTEREST:" << bby << "T.");
  Info("REGION A: " << setprecision(1) << fixed << (x1_A-x0_A)*1e+3
       << "x" << setprecision(1) << fixed << (y1_A-y0_A)*1e+3 << "mm^2." );

  // generate points at Region B
  const double x0_B = -0.047;
  const double x1_B =  0.047;
  const int    nx_B = 25;
  const double dx_B = (x1_B-x0_B) / (nx_B-1);
  const double y0_B = -0.008;
  const double y1_B =  0.008;
  const int    ny_B =  9;
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

void generate_septa_voi(VectorXd& x, VectorXd& y, VectorXd& by)
{
  // zero field region
  const double xA1 = -0.030;
  const double xA2 =  0.000;
  const int    nxA = 21;
  const double dxA = (xA2-xA1) / (nxA-1);
  const double yA1 = -0.020;
  const double yA2 =  0.020;
  const int    nyA = 21;
  const double dyA = (yA2-yA1) / (nyA-1);

  // high field region
  const double xB1 =  0.000;
  const double xB2 =  0.030;
  const int    nxB = 21;
  const double dxB = (xB2-xB1) / (nxB-1);
  const double yB1 = -0.020;
  const double yB2 =  0.020;
  const int    nyB = 21;
  const double dyB = (yB2-yB1) / (nyB-1);

  const int npoints= nxA*nyA + nxB*nyB;
  x  = VectorXd :: Zero( npoints );
  y  = VectorXd :: Zero( npoints );
  by = VectorXd :: Zero( npoints );

  int cnt = 0;

  for (int i=0; i<nxA; i++) {
    for (int j=0; j<nyA; j++) {
      x (cnt) = xA1 + dxA*i;
      y (cnt) = yA1 + dyA*j;
      by(cnt) = 2.0;
      cnt ++;
    }
  }

  for (int i=0; i<nxB; i++) {
    for (int j=0; j<nyB; j++) {
      x (cnt) = xB1 + dxB*i;
      y (cnt) = yB1 + dyB*j;
      by(cnt) = 0.0;
      cnt ++;
    }
  }
}

void problem()
{
  XFieldStraightLine* wire = new XFieldStraightLine;
  wire->CartesianCoordinate();
  wire->SetMaxOrder(10);
  wire->SetRelativePermeability(1500);
  wire->SetIronYokeRadius(0.125, -99.0);

  VectorXd x, y, by, xs, ys;
  VectorXi flag;
  Vector3d bvec;
  
  generate_wire( 1, 0.000, flag, xs, ys);
  generate_wire( 2, 0.002, flag, xs, ys);
  generate_wire( 3, 0.004, flag, xs, ys);
  generate_wire( 4, 0.006, flag, xs, ys);
  generate_wire( 5, 0.008, flag, xs, ys);
  generate_wire( 6, 0.010, flag, xs, ys);
  generate_wire( 7, 0.012, flag, xs, ys);
  generate_wire( 8, 0.014, flag, xs, ys);
  generate_wire( 9, 0.016, flag, xs, ys);
  generate_wire(10, 0.018, flag, xs, ys);
  generate_wire(11, 0.020, flag, xs, ys);
  generate_wire(12, 0.022, flag, xs, ys);
  generate_wire(13, 0.024, flag, xs, ys);
  generate_wire(14, 0.026, flag, xs, ys);
  generate_wire(15, 0.028, flag, xs, ys);
  generate_wire(16, 0.030, flag, xs, ys);
  generate_septa_voi(x, y, by);

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
  sol->SetMaxNumber(30);
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
