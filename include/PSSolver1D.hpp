#pragma once

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "PSSolver.hpp"
#include "functions.hpp"
#include <vector>
#include <array> 
#include <iostream>
#include <math.h>
#include <algorithm>


class PSSolver1D : public PSSolver {
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
     public:
          PSSolver1D();   
};
