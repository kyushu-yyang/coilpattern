#include <iostream>
#include <iomanip>
#include <fstream>
#include "XMeshHandle2.h"
#include "XSolveTSVD.h"

XSolveTSVD :: XSolveTSVD()
  : fHandle(nullptr), fNumTruncated(6),
    fFilename(""), fSave(false)
{}

XSolveTSVD :: ~XSolveTSVD()
{
  if (fHandle) delete fHandle;
}

void XSolveTSVD :: SetSaveFile(const string filename)
{
  fFilename = filename;
  fSave = true;
}

void XSolveTSVD :: Solve(const MatrixXd& A, const VectorXd& b, VectorXd& x)
{
  cout << "SOLVING INVERSE PROBLEM BY TRUNCATED SVD: AX = B" << endl;
  truncated_svd(A, b, x);

  save_each_mode(999, x);

  const double res = reconstruction(A, b, x, true);
  cout << "SUM OF RESIDES:" << setw(15) << setprecision(6) << scientific << res/b.size() << endl; 

  // save the response matrix
  save_each_mode(-999, A.row(0));
}

void XSolveTSVD :: truncated_svd(const MatrixXd& A, const VectorXd& b, VectorXd& x)
{
  /// singular value decomposition
  /// A = U*lambda*V.T
  JacobiSVD<MatrixXd> svd( A, ComputeFullU | ComputeFullV );
  VectorXd lambda = svd.singularValues();
  MatrixXd u      = svd.matrixU();
  MatrixXd v      = svd.matrixV();

  /// sum up the vector until the truncated mode
  /// x = sum( u.T*B*v/lambda )
  x   = VectorXd :: Zero(v.rows());
  VectorXd sol = VectorXd :: Zero(v.rows());

  int nmax = lambda.size();

/*
  cout << "MATRIX A:" << endl;
  cout << setprecision(5) << scientific << u << endl;
  cout << "TRANSPOSE OF MATRIX A:" << endl;
  cout << setprecision(5) << scientific << u.transpose() << endl;
  cout << "SINGULAR VALUES:" << endl;
  cout << setprecision(5) << scientific << lambda << endl;
*/

  double residual = 1000000.;
  const int num = nmax / 100;

  for (int i=0; i<nmax; i++) {
    sol = u.transpose().row(i).dot(b) * v.transpose().row(i) / lambda(i);

    //if (i<=fNumTruncated || fabs(residual/nmax)>=1e-5)
    if (i<=fNumTruncated) {
      save_each_mode(i+1, sol);
      save_each_mode(-i-1, x);
    }

    if (i==fNumTruncated) {
      reconstruction(A,b,x,true);
      save_each_mode(998, x);
    }

    residual = reconstruction(A,b,x,false);

    if (i%num==0) cout << "MODE:" << i << ", RESIDUAL:" << setprecision(5) << scientific << residual/nmax << endl;

    x += sol;
  }

  //x = svd.solve(b); 
}

double XSolveTSVD :: reconstruction(const MatrixXd& A, const VectorXd& b, const VectorXd& x, const bool print)
{
  double residual = 0.;

  VectorXd bcal = A * x;

  // print out resides
  if (print) {
    cout << "RESIDUAL OF TSVD:" << endl;
    cout << setw(6) << "POINT" << setw(15) << "B_GOAL" << setw(15) << "B_CAL" 
         << setw(15) << "B_RES" << endl; 
  }

  if (print) {
    for (int i=0; i<b.size(); i++)
      cout << setw(6) << fixed << i
           << setw(15) << setprecision(6) << scientific << b(i)
           << setw(15) << setprecision(6) << scientific << bcal(i)
           << setw(15) << setprecision(6) << scientific << b(i) - bcal(i) << endl;
  }

  VectorXd db = b-bcal;

  return db.sum();
}

void XSolveTSVD :: save_each_mode(const int mode, const VectorXd& vec)
{
  string filename = fFilename;
  filename += "_mode_";
  filename += to_string(mode);
  filename += ".vtk";

  ofstream output( filename );  
  output << "# vtk DataFile Version 2.0" << endl;
  output << "Header" << endl;
  output << "ASCII" << endl;
  output << "DATASET UNSTRUCTURED_GRID" << endl;

  // write nodes  
  output << "POINTS " << fHandle->GetNumOfNodes() << " float" << endl;
  for (int i=0; i<fHandle->GetNumOfNodes(); i++) { 
    output << setw(14) << setprecision(6) << scientific << fHandle->GetNode(i)->GetX();
    output << setw(14) << setprecision(6) << scientific << fHandle->GetNode(i)->GetY();
    output << setw(14) << setprecision(6) << scientific << fHandle->GetNode(i)->GetZ() << endl;
  }

  // write elements
  output << "CELLS " << fHandle->GetNumOfElements() << " " << fHandle->GetNumOfElements()*4 << endl;
  for (int i=0; i<fHandle->GetNumOfElements(); i++) {
    output << "3 ";
    for (int j=0; j<fHandle->GetElement(i)->GetNumOfNodes(); j++)
      output << setw(8) << fixed << fHandle->GetElement(i)->GetNode(j)->GetId();
    output << "\n";
  }

  // element type
  output << "CELL_TYPES " << fHandle->GetNumOfElements() << endl;
  for (int i=0; i<fHandle->GetNumOfElements(); i++)
    output << "5" << endl;

  // stream function
  output << "POINT_DATA " << fHandle->GetNumOfNodes() << endl;
  output << "SCALARS stream_func float" << endl;
  output << "LOOKUP_TABLE default" << endl;
  for (int i=0; i<vec.size(); i++)
    output << setw(15) << setprecision(6) << scientific << vec(i) << endl;

  output.close();
}
