#pragma once

#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <Eigen/Dense>

namespace htm {

struct PointInfo;
struct Octahedron;
struct Constraint;
struct trixel;

class HTMAsciiParser;

struct PointInfoHash
{
    size_t operator()(PointInfo const * pt) const noexcept;
};
struct PointInfoCmp
{
    bool operator()(PointInfo const *lhs, PointInfo const *rhs) const noexcept;
};

class HTM
{
 private:
  Octahedron* _octahedron;				//< base HTM Octahedron

  std::unordered_set<PointInfo*, PointInfoHash, PointInfoCmp>	_points;	//< map that reference objects by their HTMId

  std::ofstream stream;				//< Output stream to write HTM description

 public:
  /// Add a new point to the working list
  void AddPoint(const double& ra, const double& dec);

  /// Launch the creation of the HTM using the current working list
  bool CreateHTM(void);

  /// Assign a point (single operation) to the HTM
  std::string AssignPoint(PointInfo *pt);

  /// Load points from a file
  void LoadCatalog(std::string& file);

  /// Create the base octahedron
  void CreateOctahedron(void);

  /// Delete the base Octahedron
  void DeleteOctahedron(void);

  /// Get number of pairs a the radius -+ delta from the given object
  unsigned int TwoPointsCorrelation(double &radius, double &delta);

 public:
  /// Create the singleton instance if applicable and/or return a pointer to it.
  static HTM* GetInstance(void);

  /// Delete the instance of the singleton.
  static void Delete(void);

 private:

  /// Create a new trixel structure from its parent and the zone index
  // void CreateTrixelChildren(trixel *parent, unsigned int& index);


  ///
  void constraintNotInside(trixel* trixel, const Eigen::Vector3d& p, Constraint*);

  ///
  void SelectRootTrixel(PointInfo* pt);

  /// 
  inline std::pair<double, double> CalcCoordPoint(std::pair<double, double>& a, std::pair<double, double>& b);

  /// 
  inline double Scal(std::pair<double, double>& v1, std::pair<double, double>& v2) const;

  /// 
  bool CheckPointInTriangle(std::pair<double, double> A, std::pair<double, double> B, std::pair<double, double> C, std::pair<double, double> P);

  /// Check if a point is in a triangle describe by the given boundaries
  bool PointInTriangle(const double& ra, const double& dec, double* boundaries);

  /// Create a new base trixel (on of the octahedron face)
  //trixel* CreateRootTrixel(std::string HTMId);

  /// Select the first level trixel in the given octahedron
  trixel* SelectRootOctahedronTrixel(const double& ra, const double& dec);

  /// send all trixels from the given one into the given output stream
  void Display(trixel* current, std::ofstream& fstream);

  /// Free all trixels
  void FreeAllTrixels(trixel* current);

  /// BOUNDARIES
  /// Double[4] = decMin, decMax, raMin, raMax
  double* ComputeRootTrixelBounds(trixel* trixel);

  /// Compute new bounds of a trixel's child from it's parent bounds and index
  double* ComputeTrixelBounds(const double* fatherBounds, unsigned int& index, bool& reverse);

  /// get the index an point using the parent trixel bounds
  unsigned int getIndex(double* boundaries, bool& reverse, const double& ra, const double& dec);


 private:
  static HTM* _singleton; //< HTM singleton pointer.

 private:
  /// DEFAULT CTOR
  HTM(void);

  /// DEFAULT DTOR
  ~HTM(void);
};

}
