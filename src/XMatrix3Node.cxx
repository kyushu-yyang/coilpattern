#include <iostream>
#include <iomanip>
#include "XLogger.h"
#include "XMeshHandle2.h"
#include "XMatrix3Node.h"

XMatrix3Node :: XMatrix3Node()
  : fHandle(nullptr)
{
  // setup the gauss points and weights
  // n = 2
  fGaussX = MatrixXd :: Zero(4,2);
  fGaussW = VectorXd :: Zero(4);

  fGaussX(0,0) = 0.211324865; fGaussX(0,1) = 0.166666667; fGaussW(0) = 0.197168783; 
  fGaussX(1,0) = 0.211324865; fGaussX(1,1) = 0.622008467; fGaussW(1) = 0.197168783; 
  fGaussX(2,0) = 0.788675134; fGaussX(2,1) = 0.044658198; fGaussW(2) = 0.052831216; 
  fGaussX(3,0) = 0.788675134; fGaussX(3,1) = 0.166666667; fGaussW(3) = 0.052831216; 
}

XMatrix3Node :: ~XMatrix3Node()
{
  if (fHandle) delete fHandle;
}

void XMatrix3Node :: SetOrderOfGauss(const int ngauss)
{
  if (ngauss==1) {
    fGaussX = MatrixXd :: Zero(1,2);
    fGaussW = VectorXd :: Zero(1);

    fGaussX(0,0) = 1./3.; fGaussX(0,1) = 1./3.; fGaussW(0) = 0.5;
  }
  else if (ngauss==3) {
    fGaussX = MatrixXd :: Zero(9,2);
    fGaussW = VectorXd :: Zero(9);

    fGaussX(0,0) = 0.112701665; fGaussX(0,1) = 0.100000000; fGaussW(0) = 0.068464377; 
    fGaussX(1,0) = 0.112701665; fGaussX(1,1) = 0.443649167; fGaussW(1) = 0.109543004; 
    fGaussX(2,0) = 0.112701665; fGaussX(2,1) = 0.787298334; fGaussW(2) = 0.068464377; 
    fGaussX(3,0) = 0.500000000; fGaussX(3,1) = 0.056350832; fGaussW(3) = 0.061728395; 
    fGaussX(4,0) = 0.500000000; fGaussX(4,1) = 0.250000000; fGaussW(4) = 0.098765432; 
    fGaussX(5,0) = 0.500000000; fGaussX(5,1) = 0.443649167; fGaussW(5) = 0.061728395; 
    fGaussX(6,0) = 0.887298334; fGaussX(6,1) = 0.012701665; fGaussW(6) = 0.008696116; 
    fGaussX(7,0) = 0.887298334; fGaussX(7,1) = 0.056350832; fGaussW(7) = 0.013913785; 
    fGaussX(8,0) = 0.887298334; fGaussX(8,1) = 0.100000000; fGaussW(8) = 0.008696116; 
  }

  // debug 
  Info("CHANGED THE NUMBER OF GAUSS POINTS FOR INTEGRATION.");
  Info(   setw(14) << setprecision(6) << fixed << "XI_I"
       << setw(14) << setprecision(6) << fixed << "ETA_I"
       << setw(14) << setprecision(6) << fixed << "WEIGHT" );
  for (int i=0; fGaussW.size(); i++) 
    Info(   setw(14) << setprecision(6) << fixed << fGaussX(i,0)
         << setw(14) << setprecision(6) << fixed << fGaussX(i,1)
         << setw(14) << setprecision(6) << fixed << fGaussW(i) );
}

void XMatrix3Node :: SetNumOfTargetMesh(const int num)
{
  fPosTar = MatrixXd::Zero(num,3);
}

void XMatrix3Node :: SetPointAtTarget(const int i, const double x, const double y, const double z)
{
  if (i>fPosTar.rows()-1) {
    cerr << "Input index is out of range." << endl;
    throw out_of_range("input index is out of range.");
  }

  fPosTar(i,0) = x;
  fPosTar(i,1) = y;
  fPosTar(i,2) = z;
}

Vector3d XMatrix3Node :: calc_an_element(const int n_id, XElement3Node* element, Vector3d& rvec)
{
  const double dS = element->GetArea();
  Vector3d rs = element->GetPositionVector(1./3.,1./3.);
  Vector3d jn = element->GetCurrentBasisVector(n_id);
  const double rr = sqrt( pow(rvec(0)-rs(0),2) + pow(rvec(1)-rs(1),2) + pow(rvec(2)-rs(2),2) );
  Vector3d ai = jn / rr * dS;
  return ai;
}

Vector3d XMatrix3Node :: calc_bn_element(const int n_id, XElement3Node* element, Vector3d& rvec)
{
  Vector3d rs = element->GetPositionVector(1./3.,1./3.);
  const double rr = sqrt( pow(rvec(0)-rs(0),2) + pow(rvec(1)-rs(1),2) + pow(rvec(2)-rs(2),2) );
  Vector3d jn = element->GetLoopBasisVector(n_id);

  // calculate the current vector at boundary
  //int node_idx = element->GetBoundaryNodeIndex();
  //if (node_idx>0)
  //  jn += element->GetLoopBasisVector( element->GetNode(node_idx)->GetId() );

  Vector3d bi;
  
  bi(0) = (jn(1)*(rvec(2)-rs(2)) - jn(2)*(rvec(1)-rs(1))) / pow(rr,3) * 0.5;
  bi(1) = (jn(2)*(rvec(0)-rs(0)) - jn(0)*(rvec(2)-rs(2))) / pow(rr,3) * 0.5; 
  bi(2) = (jn(0)*(rvec(1)-rs(1)) - jn(1)*(rvec(0)-rs(0))) / pow(rr,3) * 0.5;

  return bi;
}

Vector3d XMatrix3Node :: calc_bn_element_gauss(const int n_id, XElement3Node* element, Vector3d& rvec)
{
  double   rr;
  Vector3d bi(0.,0.,0.);
  Vector3d rs;

  Vector3d jn = element->GetLoopBasisVector(n_id);
  
  // calculate the current vector at boundary
  //int node_idx = element->GetBoundaryNodeIndex();
  //if (node_idx>0)
  //  jn += element->GetLoopBasisVector( element->GetNode(node_idx)->GetId() );

  for (int i=0; i<fGaussW.size(); i++) {
    // setup gauss points
    rs = element->GetPositionVector( fGaussX(i,0), fGaussX(i,1) );

    // calculate |r-r'|
    rr = sqrt( pow(rvec(0)-rs(0),2) + pow(rvec(1)-rs(1),2) + pow(rvec(2)-rs(2),2) );

    // calculate field at local coordinate
    bi(0) += (jn(1)*(rvec(2)-rs(2)) - jn(2)*(rvec(1)-rs(1))) * fGaussW(i) / pow(rr,3); 
    bi(1) += (jn(2)*(rvec(0)-rs(0)) - jn(0)*(rvec(2)-rs(2))) * fGaussW(i) / pow(rr,3); 
    bi(2) += (jn(0)*(rvec(1)-rs(1)) - jn(1)*(rvec(0)-rs(0))) * fGaussW(i) / pow(rr,3); 
  }

  return bi;
}

void XMatrix3Node :: calc_an_vec(const int node_id, const double x, const double y, const double z, Vector3d& an)
{
  vector<XElement3Node*> elements = fHandle->GetSurroundingElements( node_id );
  Vector3d rvec(x, y, z);
  Vector3d an_i;

  an(0) = 0.;
  an(1) = 0.;
  an(2) = 0.;

  for (int i=0; i<elements.size(); i++) {
    an_i = calc_an_element( node_id, elements.at(i), rvec );  
    an(0) += an_i(0);
    an(1) += an_i(1);
    an(2) += an_i(2);
  }

  an *= 1e-7;
}

void XMatrix3Node :: calc_bn_vec(const int node_id, const double x, const double y, const double z, Vector3d& bn)
{
  vector<XElement3Node*> elements = fHandle->GetSurroundingElements( node_id );
  Vector3d rvec(x, y, z);
  Vector3d bn_i;

  bn(0) = 0.;
  bn(1) = 0.;
  bn(2) = 0.;

  for (int i=0; i<elements.size(); i++) {
    bn_i = calc_bn_element_gauss( node_id, elements.at(i), rvec );  
    //bn_i = calc_bn_element( node_id, elements.at(i), rvec );  
    bn(0) += bn_i(0);
    bn(1) += bn_i(1);
    bn(2) += bn_i(2);
  }

  bn *= 1e-7;
}

void XMatrix3Node :: GetFieldMatrix(MatrixXd& bxMat, MatrixXd& byMat, MatrixXd& bzMat)
{
  bxMat.conservativeResize( fPosTar.rows(), fHandle->GetNumOfNodes() );
  byMat.conservativeResize( fPosTar.rows(), fHandle->GetNumOfNodes() );
  bzMat.conservativeResize( fPosTar.rows(), fHandle->GetNumOfNodes() );

  Vector3d bij(0.,0.,0);
  cout << fPosTar.rows() << "x" << fHandle->GetNumOfNodes() << endl;
  int cnt{0};
  const int num_print = fPosTar.rows()*fHandle->GetNumOfNodes() / 100; 

#pragma omp parallel for private(bij) 
  for (int i=0; i<fPosTar.rows(); i++) {
    for (int j=0; j<fHandle->GetNumOfNodes(); j++) {
      if (cnt%num_print==0)
        cout << " CALCULATION OF FIELD AT (I,J) = (" << setw(6) << i << "," << setw(6) << j << ") ....... " 
             << setprecision(1) << fixed << cnt/static_cast<double>(fPosTar.rows())/fHandle->GetNumOfNodes()*100 << "%" << endl;

      calc_bn_vec( j, fPosTar(i,0), fPosTar(i,1), fPosTar(i,2), bij );
      bxMat(i, j) = bij(0);
      byMat(i, j) = bij(1);
      bzMat(i, j) = bij(2);
      cnt ++;
    }
  }
}
