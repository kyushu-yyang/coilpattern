#include <iostream>
#include <cmath>
#include <stdexcept>
#include "XElement3Node.h"

XElement3Node :: XElement3Node()
  : fId(-99), fSurfId(-99)
{
  fNode.push_back( nullptr );
  fNode.push_back( nullptr );
  fNode.push_back( nullptr );
}

XElement3Node :: ~XElement3Node()
{
  if (fNode.size()!=0) {
    for (auto eachnode : fNode) delete eachnode;
    fNode.clear();
  }
}

void XElement3Node :: GetNode(const int id, XNode* node)
{
  for (int i=0; i<fNode.size(); i++) {
    if (id == fNode.at(i)->GetId()) { node = fNode.at(i); break; }
  }
}

Vector3d XElement3Node :: GetShapeFunction(const double xi, const double eta)
{
  if (xi>1 || xi<0 || eta>1 || eta<0) {
    cerr << "Input parameter xi or eta is out of range." << endl;
    throw out_of_range("xi or eta is out of range.");
  }

  Vector3d N( xi, eta, 1-eta-xi );
  return N;
}

Vector3d XElement3Node :: GetPositionVector(const double xi, const double eta)
{
  Vector3d Nt = GetShapeFunction( xi, eta );
  Vector3d rr(0.,0.,0.);

  rr(0) = Nt(0)*fNode.at(0)->GetX() + Nt(1)*fNode.at(1)->GetX() + Nt(2)*fNode.at(2)->GetX();
  rr(1) = Nt(0)*fNode.at(0)->GetY() + Nt(1)*fNode.at(1)->GetY() + Nt(2)*fNode.at(2)->GetY();
  rr(2) = Nt(0)*fNode.at(0)->GetZ() + Nt(1)*fNode.at(1)->GetZ() + Nt(2)*fNode.at(2)->GetZ();

  return rr;
}

double XElement3Node :: GetArea()
{
  const double x1 = fNode.at(0)->GetX();
  const double x2 = fNode.at(1)->GetX();
  const double x3 = fNode.at(2)->GetX();

  const double y1 = fNode.at(0)->GetY();
  const double y2 = fNode.at(1)->GetY();
  const double y3 = fNode.at(2)->GetY();

  const double z1 = fNode.at(0)->GetZ();
  const double z2 = fNode.at(1)->GetZ();
  const double z3 = fNode.at(2)->GetZ();

  const double ii = (y1-y3)*(z2-z3) - (z1-z3)*(y2-y3);
  const double jj = (z1-z3)*(x2-x3) - (x1-x3)*(z2-z3);
  const double kk = (x1-x3)*(y2-y3) - (y1-y3)*(x2-x3);

  return sqrt( pow(ii,2) + pow(jj,2) + pow(kk,2) )*0.5;
}

MatrixXd XElement3Node :: GetCurrentBasisVector()
{
  MatrixXd jvec(3,3);
  
  const double x1 = fNode.at(0)->GetX();
  const double x2 = fNode.at(1)->GetX();
  const double x3 = fNode.at(2)->GetX();

  const double y1 = fNode.at(0)->GetY();
  const double y2 = fNode.at(1)->GetY();
  const double y3 = fNode.at(2)->GetY();

  const double z1 = fNode.at(0)->GetZ();
  const double z2 = fNode.at(1)->GetZ();
  const double z3 = fNode.at(2)->GetZ();

  const double A  = GetArea();

  // node 1
  jvec(0,0) = -(x2-x3)*0.5/A; 
  jvec(0,1) = -(y2-y3)*0.5/A; 
  jvec(0,2) = -(z2-z3)*0.5/A; 

  // node 2
  jvec(1,0) =  (x1-x3)*0.5/A; 
  jvec(1,1) =  (y1-y3)*0.5/A; 
  jvec(1,2) =  (z1-z3)*0.5/A; 

  // node 3
  jvec(2,0) = ((x2-x3) - (x1-x3))*0.5/A; 
  jvec(2,1) = ((y2-y3) - (y1-y3))*0.5/A; 
  jvec(2,2) = ((z2-z3) - (z1-z3))*0.5/A; 

  return jvec;
}

MatrixXd XElement3Node :: GetLoopBasisVector()
{
  MatrixXd jvec(3,3);
  
  const double x1 = fNode.at(0)->GetX();
  const double x2 = fNode.at(1)->GetX();
  const double x3 = fNode.at(2)->GetX();

  const double y1 = fNode.at(0)->GetY();
  const double y2 = fNode.at(1)->GetY();
  const double y3 = fNode.at(2)->GetY();

  const double z1 = fNode.at(0)->GetZ();
  const double z2 = fNode.at(1)->GetZ();
  const double z3 = fNode.at(2)->GetZ();

  // node 1
  jvec(0,0) = -(x2-x3); 
  jvec(0,1) = -(y2-y3); 
  jvec(0,2) = -(z2-z3); 

  // node 2
  jvec(1,0) = x1-x3; 
  jvec(1,1) = y1-y3; 
  jvec(1,2) = z1-z3; 

  // node 3
  jvec(2,0) = x2-x1; 
  jvec(2,1) = y2-y1; 
  jvec(2,2) = z2-z1; 

  return jvec;
}

Vector3i XElement3Node :: GetNodeId()
{
  Vector3i id;
  id(0) = fNode.at(0)->GetId();
  id(1) = fNode.at(1)->GetId();
  id(2) = fNode.at(2)->GetId();
  return id;
}

Vector3d XElement3Node :: GetCurrentBasisVector(const int node)
{
  // find index of the given node
  int index = -999;
  for (int i=0; i<fNode.size(); i++) {
    if (fNode.at(i)->GetId()==node) { index = i; break; }
  }

  // calculate current basis vector
  MatrixXd jvec = GetCurrentBasisVector();
  Vector3d vec( jvec(index,0), jvec(index,1), jvec(index,2) ); 

  return vec;
}

Vector3d XElement3Node :: GetLoopBasisVector(const int node)
{
  // find index of the given node
  int index = -999;
  for (int i=0; i<fNode.size(); i++) {
    if (fNode.at(i)->GetId()==node) { index = i; break; }
  }

  // calculate current basis vector
  MatrixXd jvec = GetLoopBasisVector();
  Vector3d vec( jvec(index,0), jvec(index,1), jvec(index,2) ); 

  return vec;
}

int XElement3Node :: GetBoundaryNodeIndex()
{
  int id = -999;

  // search for the node at boundary
  int cnt = 0;

  for (int i=0; i<fNode.size(); i++) {
    if (fNode.at(i)->Is_AtBoundary()) cnt++;
  }

  if (cnt>1) {
    for (int i=0; i<fNode.size(); i++) {
      if (!fNode.at(i)->Is_AtBoundary()) id = i;
    }
  }

  return id;
}
