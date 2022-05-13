#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include "XVtkFile.h"
#include "XMeshHandle2.h"
#include "XElement3Node.h"
#include "XMatrix3Node.h"
#include "XSolveTSVD.h"
#include "XLogger.h"

using namespace std;
using namespace Eigen;


void write_Amat(const char* filename, const MatrixXd& mat)
{
  ofstream input( filename );

  for (int i=0; i<mat.rows(); i++) {
    for (int j=0; j<mat.cols(); j++) 
      input << setw(15) << setprecision(6) << scientific << mat(i,j);
    input << "\n";
  }

  input.close();
}

void target_point(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& by)
{
  const double rho0   = 1.89;
  const double phi0   = M_PI/32-M_PI/80;
  const double phi1   = M_PI/32+M_PI/80;
  const int    nphi   = 10;
  const double theta0 = 0.;
  const double theta1 = 2*M_PI;
  const int    ntheta = 10;
  const double r0     = 0.004;
  const double r1     = 0.030;
  const int    nr     = 2;
  const double bnorm  = 2.4;

  double xx, yy, zz, r, theta, phi;

  x.clear();
  y.clear();
  z.clear();

  for (int i=0; i<nphi; i++) {
    phi = phi0 + (phi1-phi0)*i/(nphi-1);
    for (int j=0; j<nr; j++) {
      r = r0 + (r1-r0)*j/(nr-1);
      for (int k=0; k<ntheta; k++) {
        theta = theta0 + (theta1-theta0)*k/(ntheta-1);
        xx = (rho0 + r*cos(theta)) * cos(phi);
        yy = r*sin(theta);
        zz = -(rho0 + r*cos(theta)) * sin(phi);
 
        x .push_back(xx);
        y .push_back(yy);
        z .push_back(zz);
        by.push_back(bnorm);
      }
    }
  }
}

void dsv_point(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& bz)
{
  const double phi0   = 0.0;
  const double phi1   = 2*M_PI;
  const int    nphi   = 8;
  const double theta0 = 0.;
  const double theta1 = M_PI;
  const int    ntheta = 8;
  const double r0     = 0.001;
  const double r1     = 0.180;
  const int    nr     = 5;
  const double z0     = 1.0;
  const double grad   = 1.0;

  double xx, yy, zz, r, theta, phi;

  x.clear();
  y.clear();
  z.clear();

  x .push_back(0.);
  y .push_back(0.);
  z .push_back(z0);
  bz.push_back(0.);

  for (int i=0; i<nr; i++) {
    r = r0 + i*(r1-r0)/(nr-1);
    for (int j=0; j<ntheta; j++) {
      theta = theta0 + j*(theta1-theta0)/(ntheta-1);
      for (int k=0; k<nphi; k++) {
        phi = phi0 + (phi1-phi0)*k/(nphi-1);

        xx = r * sin(theta) * cos(phi);
        yy = r * sin(theta) * sin(phi);
        zz = r * cos(theta) + z0;

        x .push_back(xx);
        y .push_back(yy);
        z .push_back(zz);
        bz.push_back(grad*xx);
      }
    }
  }
}

void test(const char* filename)
{
  XVtkFile* file = new XVtkFile( filename );
  //file->PrintNode();
  //file->PrintElement();

  XMeshHandle2* handle = new XMeshHandle2;
  handle->SetVtkFile( file );
  //handle->Print();
  delete file;

  //XMeshHandle* handle = new XMeshHandle( "../input/curvedGeo.vtk" );
  const int numOfNodes    = handle->GetNumOfNodes();  
  const int numOfElements = handle->GetNumOfElements();  

  cout << "NUMBER OF NODES   :" << numOfNodes    << endl;
  cout << "NUMBER OF ELEMENTS:" << numOfElements << endl;

  /**********************************************************/
  XMatrix3Node* mat = new XMatrix3Node;
  mat->SetMeshInfo(handle);

/*
  const double rho0   = 1.89;
  const double phi0   = 0.;
  const double phi1   = M_PI/4;
  const int    nphi   = 3;
  const double theta0 = 0.;
  const double theta1 = 2*M_PI;
  const int    ntheta = 25;
  const double r0     = 0.001;
  const double r1     = 0.030;
  const int    nr     = 10;

  vector<double> xt;
  vector<double> yt;
  vector<double> zt;
  double r, phi, theta;

  for (int i=0; i<nphi; i++) {
    phi = phi0 + (phi1-phi0)*i/(nphi-1);
    for (int j=0; j<ntheta; j++) {
      theta = theta0 + (theta1-theta0)*j/(ntheta-1);
      for (int k=0; k<nr; k++) {
        r = r0 + (r1-r0)*k/(nr-1);

        xt.push_back( (rho0+r*cos(theta)) * cos(phi) );
        yt.push_back( r*sin(theta) );
        zt.push_back(-(rho0+r*cos(theta)) * sin(phi) );
      }
    }
  }
*/

  vector<double> xt;
  vector<double> yt;
  vector<double> zt;
  vector<double> byy;

  target_point(xt, yt, zt, byy);
  //dsv_point(xt, yt, zt, byy);
 
  mat->SetNumOfTargetMesh( xt.size() );

  VectorXd b = VectorXd :: Ones( byy.size() );

  for (int i=0; i<zt.size(); i++) {
    mat->SetPointAtTarget(i, xt.at(i), yt.at(i), zt.at(i) ); 
    b(i) = byy.at(i);
    //cout << i << " " << xt.at(i) << " " << yt.at(i) << " " << zt.at(i) << endl;
  }

  cout << "******** CALCULATION OF DENSE MATRIX *********" << endl;

  MatrixXd bx_mat, by_mat, bz_mat;
  mat->GetFieldMatrix( bx_mat, by_mat, bz_mat );

  write_Amat( "matrixBx.txt", bx_mat );
  write_Amat( "matrixBy.txt", by_mat );
  write_Amat( "matrixBz.txt", bz_mat );
  cout << " DENSE MATRIX IS READY" << endl;
  /**********************************************************/
  VectorXd x = VectorXd :: Zero( by_mat.cols() ); 

  cout << " SOLVING THE PROBLEM: AX = B" << endl;

  XSolveTSVD* solver = new XSolveTSVD;
  solver->SetMeshHandle(handle);
  solver->SetTruncatedNumber(120);
  solver->SetSaveFile("res_tsvd");
  solver->Solve( by_mat, b, x );
  //solver->Solve( bz_mat, b, x );

  /*
  //VectorXd xx = by_mat.colPivHouseholderQr().solve(b);
  JacobiSVD<MatrixXd> svd( by_mat, ComputeFullU | ComputeFullV );
  VectorXd xx = svd.solve(b);
  VectorXd lam= svd.singularValues();
  MatrixXd uu = svd.matrixU();
  MatrixXd vv = svd.matrixV();

  cout << setprecision(6) << scientific << lam << endl;
  cout << setprecision(6) << scientific << uu  << endl;

  cout << setprecision(6) << scientific << uu.transpose().row(2).dot(b)  << endl;
  VectorXd sol = uu.transpose().row(2).dot(b) * vv.row(2) / lam(2);

  //cout << xx << endl;
  //cout << by_mat * xx << endl;
  */
}

int main(int argc, char** argv)
{
  XLogger* log = XLogger::GetInstance();
  log->Start(XLogger::DEBUG, "debug_.log");

  try {
    test(argv[1]);
  }
  catch (invalid_argument& except) {
    cerr << except.what() << endl;
    return EXIT_FAILURE;
  }

  log->Stop();

  return EXIT_SUCCESS;
}
