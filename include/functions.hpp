#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED
 
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <vector>
#include <array> 
#include <iostream>

class Vec3d;

/**
 * This class is used to calculate the surface of polygons
 */
class Surface {
private:
  std::vector<Vec3d> vertex;
  int numbers;

public:
  void init(int length);
  ///This function stores the vertices
  void import_vertex(std::vector<Vec3d> data);
  ///This function calculates the surface of the given polygon
  double get_surface(Vec3d n);
};

/**
 * This class implements a 2d vector
 */
class Vec2d{
    private: 
    
    public:

        std::array<double, 2> koord;
        ///This function multiplies the vector with scalar
        Vec2d operator*(double scalar);
        ///This function add another vector
        Vec2d operator+(Vec2d p2);
        ///This function subtracts another vector
        Vec2d operator-(Vec2d p2);
        void operator=(std::array<double, 2> p2);
        void operator=(std::vector<double> p2);
        ///This function calculates the scalar product
        double operator*(Vec2d p2); 
    
};

/**
 * This class implements a 3d vector
 */
class Vec3d{
    private: 
    
    public:
    
    std::array<double, 3> koord;
    
    ///This function multiplies the vector with scalar
    Vec3d operator*(double scalar);
    ///This function add another vector
    Vec3d operator+(Vec3d p2);
    ///This function subtracts another vector
    Vec3d operator-(Vec3d p2);
    void operator=(std::array<double, 3> p2);
    void operator=(std::vector<double> p2);
    ///This function calculates the scalar product
    double operator*(Vec3d p2);
    ///This function calculates the norm of the vector
    double norm();
    ///This function calculates the cross product
    Vec3d operator^(Vec3d B);
};

/**
 * This class implements a solver for a 2x2 matrix
 */
class Matrix2d{
    private:
        std::array<double, 4> koeff = {};  //a11, a12, a21, a22
    public:
        ///This function sets one coefficient of the matrix
        void setkoeff(long row, long column, double value);
        ///This function solves the system of linear equations
        std::array<double,2> solve(std::array<double, 2> vec);
};

/**
 * This class implements a solver for a 3x3 matrix
 */
class Matrix3d{
    private:
        std::array<double, 9> koeff = {};  //a11, a12, a13, a21, a22, a23, a31, a32, a33
    public:
        ///This function sets one coefficient of the matrix
        void setkoeff(long row, long column, double value);
         ///This function solves the system of linear equations
        std::array<double,3> solve(std::array<double, 3> vec);
};


/**
 * This class stores a delaunay meh on a plane and calculates the surface of each control volume
 */
class mesh2d{
    public:
        ///stores the indices of the points for each triangle
        std::vector<std::array<long, 3> > triangles = {};
        ///stores the circumcenter for each triangle
        std::vector<std::array<double, 3> > circumcenter = {};
        ///stores the points on the plane
        std::vector<std::array<double,3> > points = {};
        ///stores the points indices
        std::vector<long> points_indices = {};
        ///one point which lies on the plane
        Vec3d surfacepoint;
        ///outward pointing normalvector
        Vec3d normalvector;
        ///This function add new triangles and points to the plane
        void init(std::array<long, 3> tr, std::array<double,3> p1,
                    std::array<double,3> p2, std::array<double,3> p3);
        ///This function calculates the control volume of one point
        double getvolume(long center);
};

/**
 * This struct is needed to store one point and the index of the point before and after.  Only 2D points are stored
*/

struct List2D{
    long neighbor_prev = -1, neighbor_next = -1;
    std::array<double,2> p = {};
    ///This function stores value which is from type Vec2d
    void operator=(Vec2d value) {
        p[0] = value.koord[0];
        p[1] = value.koord[1];
    }
    ///This function stores value which is from type std::vector<double>
    void operator=(std::vector<double> value) {
        p[0] = value[0];
        p[1] = value[1];
    }
};

/**
 * This struct is needed to store one point and the index of the point before and after.  Only 3D points are stored
*/

struct List3D{
    long neighbor_prev = -1, neighbor_next = -1;
    std::array<double,3> p = {};
    ///This function stores value which is from type Vec3d
    void operator=(Vec3d value) {
        p[0] = value.koord[0];
        p[1] = value.koord[1];
        p[2] = value.koord[2];

    }
    ///This function stores value which is from type std::vector<double>
    void operator=(std::vector<double> value) {
        p[0] = value[0];
        p[1] = value[1];
        p[2] = value[2];

    }
};


#endif // FUNCTIONS_H_INCLUDED
