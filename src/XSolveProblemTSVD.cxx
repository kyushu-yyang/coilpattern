#include <iostream>
#include <iomanip>
#include "XLogger.h"
#include "XSolveProblemTSVD.h"

XSolveProblemTSVD :: XSolveProblemTSVD()
  : fNmax(10)
{}

void XSolveProblemTSVD :: SetVolumeOfInterest(const VectorXd& x, const VectorXd& y, const VectorXd& b) 
{
  fPos = MatrixXd :: Zero( x.size(), 2 );
  fBtar= VectorXd :: Zero( x.size() );

  for (int i=0; i<x.size(); i++) {
    fPos(i,0) = x(i);
    fPos(i,1) = y(i);
    fBtar(i)  = b(i);
  }
}

void XSolveProblemTSVD :: SetSourcePosition(const VectorXi flag, const VectorXd& xs, const VectorXd& ys)
{
  fFlag   = flag;
  fPosSrc = MatrixXd :: Zero( xs.size(), 2 );

  for (int i=0; i<xs.size(); i++) {
    fPosSrc(i,0) = xs(i);
    fPosSrc(i,1) = ys(i);
  }
}

void XSolveProblemTSVD :: Solve(const MatrixXd& A)
{
  VectorXd x;

  // solve the problem
  truncated_SVD(A, fBtar, x);

  // reconstruction
  cout << "RECONSTRUCTION OF FIELD" << endl;
  cout << setw(14) << "X"
       << setw(14) << "Y"
       << setw(14) << "BTAR[T]"
       << setw(14) << "BREC[T]"
       << setw(14) << "BRES[T]" << endl;

  VectorXd brec = A * x;

  for (int i=0; i<fBtar.size(); i++)
    cout << setw(14) << setprecision(6) << scientific << fPos(i,0)
         << setw(14) << setprecision(6) << scientific << fPos(i,1)
         << setw(14) << setprecision(6) << scientific << fBtar(i)
         << setw(14) << setprecision(6) << scientific << brec(i)
         << setw(14) << setprecision(6) << scientific << fBtar(i) - brec(i) << endl;

  // result of ideal current distribution
  cout << "IDEAL DISTRIBUTION OF COIL CURRENT." << endl;
  for (int i=0; i<x.size(); i++)
    cout << setw(15) << setprecision(6) << scientific << x(i) << endl;
}

void XSolveProblemTSVD :: Write(const string filename)
{
  ofstream output( filename );

  output << "MAXIMUM MODE NUMBER:" << setw(4) << fNmax << endl;

  output << "MODE STRENGTH " << fPtg.size() << endl;
  output << setw(6) << "MODE" << setw(14) << "STRENGTH" << endl;
  for (int i=0; i<fPtg.size(); i++)
    output << setw(6) << fixed << i
           << setw(14) << setprecision(6) << scientific << fPtg(i) << endl; 

  output << "MODE OF MAGNETIC FIELD " << fBtar.size() << endl;
  output << setw(14) << "X" << setw(14) << "Y";
  for (int i=0; i<fNmax; i++)
    output << setw(15) << i;
  output << "\n";

  for (int i=0; i<fBtar.size(); i++) {
    output << setw(14) << setprecision(6) << scientific << fPos(i,0);
    output << setw(14) << setprecision(6) << scientific << fPos(i,1);

    for (int j=0; j<fNmax; j++)
      output << setw(15) << setprecision(6) << scientific << fBrec(j,i);

    output << "\n";
  }

  output << "MODE OF CURRENT " << fIvec.rows() << " " << fIvec.cols() << endl;

  if (fFlag.size()!=0)
    output << setw(5 ) << "FLAG"
           << setw(15) << "XS"
           << setw(15) << "YS";

  for (int i=0; i<fNmax; i++)
    output << setw(15) << i;
  output << "\n";

  for (int i=0; i<fIvec.cols(); i++) {
    if (fFlag.size()!=0)
      output << setw( 5) << fixed << fFlag(i)
             << setw(15) << setprecision(6) << scientific << fPosSrc(i,0)
             << setw(15) << setprecision(6) << scientific << fPosSrc(i,1);

    for (int j=0; j<fNmax; j++)
      output << setw(15) << setprecision(6) << scientific << fIvec(j,i);

    output << "\n";
  }

  output.close();
}

void XSolveProblemTSVD :: truncated_SVD(const MatrixXd& A, const VectorXd& b, VectorXd& x)
{
  // singular value decomposition
  // A = U*lambda*V.T
  JacobiSVD<MatrixXd> svd( A, ComputeFullU | ComputeFullV );
  VectorXd lambda = svd.singularValues();
  MatrixXd u      = svd.matrixU();
  MatrixXd v      = svd.matrixV();

  // check the maximum truncated number
  if (fNmax > b.size()) fNmax = b.size();

  // sum up the vector until the truncated mode
  // x = sum( u.T*B*v/lambda )
  x = VectorXd :: Zero(v.rows());
  fIvec = MatrixXd :: Zero( fNmax, x.size() );
  fBrec = MatrixXd :: Zero( fNmax, b.size() );
  fPtg  = VectorXd :: Zero( fNmax );
  VectorXd I, brec;

  double res, P_TG;
  const int num = fNmax / 100;
  const int npts= b.size(); 

  // sum the each mode upto the truncated number 
  for (int i=0; i<fNmax; i++) {
    P_TG = u.transpose().row(i).dot(b) / sqrt(npts);
    I    = sqrt(npts) * v.transpose().row(i) * P_TG / lambda(i);
    brec = sqrt(npts) * P_TG * u.col(i); 

    Info( "MODE-" << i << ", P_TG:" << setw(14) << setprecision(6) << scientific << P_TG << "T." );

    // fill the matrix
    for (int j=0; j<fIvec.cols(); j++) fIvec(i,j) = I(j);
    for (int j=0; j<fBrec.cols(); j++) fBrec(i,j) = brec(j);
    fPtg(i) = P_TG;

    x += I;
  }

}
