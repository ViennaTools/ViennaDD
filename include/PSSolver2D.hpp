#ifndef DIFFUSION2D_H_INCLUDED
#define DIFFUSION2D_H_INCLUDED

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "PSSolver.hpp"
#include "functions.hpp"
#include <vector>
#include <array> 
#include <iostream>
#include <math.h>
#include <algorithm>

/**
 * This class is inherited from the class PSSolver and calculates the control volumes in 2D, surfaces and the circumcenters in 2D
 */
class PSSolver2D : public PSSolver {
private:
    
     ///This function calculates the circumcenter for every triangle in the delaunay mesh
     virtual void create_mesh();
     ///This function calculates the control surface of the point indices center and neighbor
     virtual double getsurface(long center, long neighbor);
     ///This function calculates the control volume in 2D
     virtual double getvolume(long center);
     ///This function calculates the distance between the points center and neighbor
     virtual double getdistance(long center, long neighbor);

     virtual void scalepoints(std::vector<std::vector<double>> &grid, double factor);

     ///number of points in x direction
     long Nx = 0;
     ///number of points in y direction
     long Ny = 0;
     
     public:
          PSSolver2D();
          ///This function imports an rectangular mesh and save it as a triangular mesh
          void import_rectangular_mesh(std::vector<std::array<long,4> > mesh, long Nx, long Ny);
          ///This function sets the number of points in each direction, used for rectangular meshes
          void set_number_of_points(long Nx, long Ny);
          ///This function calculates the volume of the device and prints it out
          void test();
          ///This function solves the poisson equation for an rectangular mesh
          std::vector<double> poisson_rect();
          
};

#endif // DIFFUSION1D_H_INCLUDED
