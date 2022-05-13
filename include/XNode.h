#ifndef XNODE_HH
#define XNODE_HH

#include <Eigen/Dense>

using namespace Eigen;

/// @file   XNode.h
/// @author Ye Yang (QST)
/// @date   05.13.2022

class XNode
{
  public:
    XNode() { 
      fXYZ(0) = 0.;
      fXYZ(1) = 0.;
      fXYZ(2) = 0.; 
    }

    ~XNode() {}

    void SetId(const int num) { fId = num; }
    int  GetId() const { return fId; }

    void SetPoint(const double xx, const double yy, const double zz) { 
      fXYZ(0) = xx;
      fXYZ(1) = yy;
      fXYZ(2) = zz;
    }
    Vector3d GetPoint() const { return fXYZ; }

    void   SetX(const double xx) { fXYZ(0)=xx; }
    void   SetY(const double yy) { fXYZ(1)=yy; }
    void   SetZ(const double zz) { fXYZ(2)=zz; }

    double GetX() const { return fXYZ(0); }
    double GetY() const { return fXYZ(1); } 
    double GetZ() const { return fXYZ(2); }

    double GetTheta () { return atan2( fXYZ(1), fXYZ(0) ); }
    double GetRadius() { return sqrt(pow(fXYZ(0),2) + pow(fXYZ(1),2)); }

    Vector3d GetPositionVector(const double x, const double y, const double z) {
      Vector3d vec(0.,0.,0.);
      vec(0) = x - fXYZ(0);
      vec(1) = y - fXYZ(1);
      vec(2) = z - fXYZ(2);
      return vec;
    }
    double GetDistance(const double x, const double y, const double z) {
      Vector3d vec = GetPositionVector(x,y,z);
      return vec.norm();
    }

    void   AtBoundary(const int id) {
      fNodeAtBC = true; 
      fBcId     = id; 
    }
    bool   Is_AtBoundary() const { return fNodeAtBC; }
    int    GetBoundaryId() const { return fBcId; }

  private:
    int      fId{-99};
    Vector3d fXYZ;
    int      fBcId{-99};
    bool     fNodeAtBC{false};
};

#endif
