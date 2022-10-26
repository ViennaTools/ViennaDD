#include "PSSolver1D.hpp"


void PSSolver1D::create_mesh() {
    ut = kb*T/q;
    ni = sqrt(Nc*Nv*exp((-(Ec-Ev)/ut)));
}

double PSSolver1D::getsurface(long center, long neighbor) {
    if(center > 0 && neighbor > 0) return 1;
    else return 1;
}

double PSSolver1D::getvolume(long center) {
    double a = points[center+1][0];
    double b = points[center-1][0];
    return std::abs(a-b)/2;
}

double PSSolver1D::getdistance(long center, long neighbor) {
    double a = points[center][0];
    double b = points[neighbor][0];
    return std::abs(a-b);
}
          
     
PSSolver1D::PSSolver1D() {
    dim = 2;
}


void PSSolver1D::scalepoints(std::vector<std::vector<double>> &grid, double factor) {
    for(auto it = grid.begin(); it !=grid.end(); it++) {
        for(auto it2 = it->begin(); it2 != it->end(); it2++) {
            *it2 = *it2*factor;
        }
    }
}
