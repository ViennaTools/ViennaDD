#ifndef DIFFUSION3D_H_INCLUDED
#define DIFFUSION3D_H_INCLUDED


#include <iostream>
#include "Eigen/IterativeLinearSolvers"
#include <math.h>
#include <fstream>
#include <algorithm>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "PSSolver.hpp"
#include <vector>
#include <array> 
#include <map>
#include "functions.hpp"


/**
 * This class is inherited from the class PSSolver and calculates the control volumes in 3D, surfaces and the circumcenters in 3D
 */
class PSSolver3D : public PSSolver{
private:
  
  ///saves the boundary surfaces and their normalvector
  std::vector<mesh2d> simple_boundaries;
  ///initialize simple_boundaries
  void init_surface();
   ///This function calculates the circumcenter for every triangle in the delaunay mesh
  virtual void create_mesh();

  ///all boundary points  
  std::vector<std::array<double,3> > polygonpoints = {};
  ///saves the inces of the polygonpoints which coresspond to a boundary plane
  std::vector<std::vector<long> > surfacepolygons = {};
  ///One point which is somewhere in the middle of the simulation area
  Vec3d middlepoint;

  ///This function calculates the control surface of the point indices center and neighbor
  virtual double getsurface(long center, long neighbor1);
  ///This function calculates the control volume in 3D
  virtual double getvolume(long center);
  ///This function calculates the distance between the points center and neighbor
  virtual double getdistance(long center, long neighbor);
  virtual void scalepoints(std::vector<std::vector<double>> &grid, double factor);

  ///constists of the indices of surfacepolygons for each simulation. 
  std::map<const std::string, std::vector<long> > simple_segments;

  public:
    ///saves points in polygonpoints
    void set_polygonpoints(std::vector<std::array<double,3> > points);
    ///saves points in surfacepolygon
    void set_surfacepolygon(std::vector<std::vector<long> >  points);
    ///add boundary and name to simple_segments
    void add_simple_boundary(std::vector<long> boundary, std::string name);
    PSSolver3D();
    ///This function calculates the volume of the device and prints it out
    void test();
};



#endif // DIFFUSION3D_H_INCLUDED
