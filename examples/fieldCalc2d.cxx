#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <Eigen/Dense>
#include "XFieldStraightLine.h"
#include "XLogger.h"

using namespace std;
using namespace Eigen; 

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
    double xx, yy, phi;

    if (!(iss >> flag >> phi >> xx >> yy >> curr))
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
  // generate points at Region A
  const double x0_A = -0.130;
  const double x1_A =  0.130;
  const int    nx_A = 60;
  const double dx_A = (x1_A-x0_A) / (nx_A-1);
  const double y0_A = -0.130;
  const double y1_A =  0.130;
  const int    ny_A = 60;
  const double dy_A = (y1_A-y0_A) / (ny_A-1);
  const double bby  = 0.0;

  Info("FIELD AT VOLUME OF INTEREST:" << bby << "T.");
  Info("REGION A: " << setprecision(1) << fixed << (x1_A-x0_A)*1e+3
       << "x" << setprecision(1) << fixed << (y1_A-y0_A)*1e+3 << "mm^2." );

  /*
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
  */

  // generate points
  //int npoints = nx_A*ny_A + nx_B*ny_B;
  int npoints = nx_A*ny_A;
  x = VectorXd :: Zero( npoints );
  y = VectorXd :: Zero( npoints );
  bx= VectorXd :: Zero( npoints );
  by= VectorXd :: Zero( npoints );
  double xx, yy;
  int cnt = 0;

  for (int i=0; i<nx_A; i++) {
    xx = x0_A + dx_A*i;
    for (int j=0; j<ny_A; j++) {
      yy = y0_A + dy_A*j;
      x (cnt) = xx;
      y (cnt) = yy;
      bx(cnt) = bby;
      by(cnt) = bby;
      cnt ++;
    }
  }

  /*
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
  */
}

void problem(const string filename)
{
  XFieldStraightLine* wire = new XFieldStraightLine;
  wire->CartesianCoordinate();
  wire->SetMaxOrder(6);
  wire->SetRelativePermeability(1500);
  wire->SetIronYokeRadius(0.125, -99.0);

  Vector3d bvec;
  VectorXd x, y, by, bx;
  vector<double> xs, ys, Is;
  const double Inorm = 180.;

  load_file( filename, xs, ys, Is ); 
  
  generate_voi(x, y, bx, by);

  for (int i=0; i<x.size(); i++) {
    for (int j=0; j<xs.size(); j++) {
      wire->SetPoint( xs.at(j), ys.at(j), Is.at(j)*Inorm );
      bvec = wire->GetMagneticField( x(i), y(i) );
      bx(i) += bvec(0);
      by(i) += bvec(1);
    }
  }

  cout << setw(15) << "X[m]" << setw(15) << "Y[m]" << setw(15) << "BY[T]" << endl;

  for (int i=0; i<x.size(); i++) 
    cout << setw(15) << setprecision(6) << scientific << x (i) 
         << setw(15) << setprecision(6) << scientific << y (i) 
         << setw(15) << setprecision(6) << scientific << bx(i) 
         << setw(15) << setprecision(6) << scientific << by(i) << endl; 
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
