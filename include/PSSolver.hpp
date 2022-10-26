#pragma once 

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <vector>
#include <array> 
#include <string>
#include <iostream>
#include <map>
#include "functions.hpp"
#include <math.h>
#include <algorithm>
#include <iterator>
#include "Eigen/IterativeLinearSolvers"

/**
 * This class solves the Drift-Diffusion-Modell in 1D
 */
class PSSolver {
    protected:

    double mu_electron = 1200, mu_holes = 400, epsilon_r = 11.9, ut = 0.025;
    double  ni = 1e10;
    ///in eV
    double Ec = 5.6198749999999997e-01;
    ///in eV
    double Ev = -5.6198749999999997e-01;
    double T = 300, Nc =  2.8541782947081875e+19,
          Nv = 2.6998702134944203e+19, kb = 1.38064852e-23, q = 1.602176634e-19;

    double Ei;
    double intrinsic_debye;
    double steps = 1;
    double Et = 1.1999999999999997e+00;
    double n1,p1;
    double elec_mobility = 500, hole_mobility = 50;
    double rspond = 1e-10;
    double Cn = 1.1e-30;
    double Cp = 3e-31;
    bool solve_same_time = true;
    bool SHR = false;
    bool Auger = false;
    bool Spond = false;
    double elec_lifetime = 1.0779393682159489e-10, hole_lifetime =  1.4359058527257259e-11;

    Eigen::VectorXd pot;
    std::vector<double> Efn;
    std::vector<double> Efp;

    long dim;
    ///maximum number of gummel iterations
    long max_gummel_iterations = 100;
    ///maximum number of newton iterations
    long max_newton_iterations = 100;
    double epsilongummel = 1e-8, epsilonnewton = 1e-12;
    ///stores the current electron concentration
    Eigen::VectorXd elec;
    ///stores the current hole concentration
    Eigen::VectorXd holes;  

    ///bernouli function for the Scharfetter-Gummel scheme
    double bernoulifunction(double Vi, double Vj);
    ///This function calculates the hole concentration with the Boltzman statistics
    double holes_fermi(double voltage);
    ///This function calculates the electron concentration with the Boltzman statistics
    double elec_fermi(double voltage);
    
    ///This functions solves the poisson equation which is optained from a fixedpoint iteration
    Eigen::VectorXd fixpoint(Eigen::VectorXd oldpot);
    ///This function performs one newtoniteration for the nonlinear poisson equation in equillibrium
    Eigen::VectorXd newtoniteration(Eigen::VectorXd potold);
    ///This function calculates the holeconcentration with the Scharfetter-Gummel 
    Eigen::VectorXd holeconcentration(Eigen::VectorXd potential);
    ///This function calculates the electronconcentration with the Scharfetter-Gummel scheme
    Eigen::VectorXd electronconcentration(Eigen::VectorXd potential);
    ///This function returns the max. error of the nonlinear poisson equation
    double error(Eigen::VectorXd potold);   

    ///stores the indices of the triangles which are part of a simulation segment
    std::map<const std::string, std::vector<long> > semiconductor_segments;
    ///stores the indices of the triangles which are part of a contact segment
    std::map<const std::string, std::vector<long> > contact_segments;
    ///stores the indices of the points which are part of a contact segment
    std::map<const std::string, std::vector<long> > contact_segments_points;
    
    ///saves the dirichlet boundary conditions
    std::map<const long, double> dirichlet;
    ///saves the neumann boundary conditions
    std::map<const long, double> neumann;
    
    ///stores the control volume for each point in the simulation grid
    std::vector<double> volumes;
    ///stores the control surfaces
    std::map<long, std::map<long,double>> surfaces;
    ///stores the indices of the points which are used in the simulation grid
    std::vector <long> indices;
    ///stores all points
    std::vector<std::vector<double>> points = {};
    ///stores the indices of the points which correspond to the delaunay mesh
    std::vector<std::vector<long> > triangles = {};
    ///stores the circumcenter for each triangle 
    std::vector<std::vector<double> > circumcenter = {}; 

    ///stores the dopand concentration for each point
    std::vector<double> dopand;
    ///stores the aceptor concentration for each point
    std::vector<double> aceptor; 
    ///stores the charge density for each point
    std::vector<double> charge_density;
    /// This function returns all neighbors for the center node
    std::vector<long> find_neighbors(long center);

    virtual void scalepoints(std::vector<std::vector<double>> &grid, double factor);
    virtual void create_mesh();
    virtual double getsurface(long center, long neighbor);
    virtual double getvolume(long center);
    virtual double getdistance(long center, long neighbor);
    /// This function solves the conservation law of electrons for the quasi-fermi level
    void electronconcentrationfermi();
    /// This function solves the conservatin law of holes for the quasi-fermi level
    void holeconcentrationfermi();

    ///This function solves the scaled nonlinear poisson equation
    Eigen::VectorXd scaled_nonlinear_poisson(Eigen::VectorXd potold);

public:
    ///This function solves the DD-model in thermodynamic equillibrium
    void calculateBuiltIn();
    /// This function solves the DD-model by solving the continuity equations for the quasi-fermi levels
    void solveDDQuasiFermi();
    ////// This function solves the DD-model by solving the continuity equations directly
    void solveDDClassic(); 
    ///This function sets the stepwide  
    void setSteps(double _steps);
    ///This function solves the poisson equation 
    std::vector<double> solvePoisson();
    ///This function sets the break criteria for the gummel method
    void setEpsilonGummel(double eps);
    ///This function sets the relativ permitivity
    void setPermitivity(double perm);
    ///This function sets the Temperature in Kelvin
    void setTemperatur(double temp);
    ///This function sets the electron mobility
    void setElectronMobility(double mob);
    ///This function sets the hole mobility
    void setHoleMobility(double mob);
    /// This function sets the factor of the sponatanous recombination
    void setSpontanRecomb(double _rspond);
    /// This function sets the Cn coefficient for the auger recombination
    void setCn(double _Cn);
    /// This functin sets the Cp coefficient of the auger recombination
    void setCp(double _Cp);
    ///This function sets the break criteria for the newton method
    void setEpsilonNewton(double eps);
    ///This function sets the electron lifetime
    void setElecLifetime(double time);
    ///This function sets the hole lifetime
    void setHoleLifetime(double time);
    ///This function sets the Trap Energy
    void setTrapEnergy(double energy);
    ///This function activates the SHR Recombination in the continuity equations
    void setSHR(bool status);
    ///This function activates the Auger Recombination in the continuity equations
    void setAuger(bool status);
    ///This function activates the spontanous Recombination in the continuity equations
    void setSpond(bool status);
    ///This function sets Nc in 1/(cm³)
    void setElecDensityofStates(double value);
    ///This function sets Nv in 1/(cm³)
    void setHolesDensityofStates(double value);
    ///This functions sets the charge density
    void setChargeDensity(std::vector<double> charge);
    /// This function sets the energy of the conduction band
    void setConductionBandEnergy(double _Ec);
     /// This function sets the energy of the valence band
    void setValenceBandEnergy(double _Ev);
    ///This function adds a contact segment to the simulation grid
    void addContactSegment(std::vector<long> segment, std::string name);
    ///This function adds points to a contact segment
    void addContactPoints(std::vector<long> points, std::string name);
    ///This function add an semiconductor segment to the simulation grid
    void addSemiconductorSegment(std::vector<long> data, std::string name);
    ///This function contact all points to an Delaunay mesh
    void setCellGrid(std::vector<std::vector<long> > data); 
    ///This function sets the points for the simulation Grid
    void setSimulationGrid(std::vector<std::vector<double> > points);
    ///This function sets the donor concentration on each point
    void setDonorConcentration(std::vector<double>  concentration);
     ///This function sets the aceptor concentration on each point
    void setAcceptorConcentration(std::vector<double>  concentration);
    ///This function sets the maximum number of gummel iterations
    void setNumberGummelIterations(long num);
    ///This function sets the maximum number of newton iterations
    void setNumberNewtonIterations(long num);
    ///This function applies an voltage to an existing contact segment
    void applyVoltageSegment(double voltage, std::string name);
    ///This function deletes the applied voltage on a contact segment
    void deleteVoltageSegment(std::string name);
    ///This function applies voltage to a single point
    void applyVoltagePoint(double voltage, long index);
    ///This function deletes the applied voltage of a single point
    void deleteVoltagePoint(long index);
    ///This function add nemann boundary conditions to an existing segment 
    void addNeumannSegment(double value, std::string name);
    ///This function deletes neumann boundary conditions of an existing segment
    void deleteNeumannSegment(std::string name);
    ///This function adds neumann boundary conditions to an existing point
    void addNeumannPoint(double value, long index);
    ///This function deletes neumann bounday condition of an existing point
    void deleteNeumannPoint(long index);
    /// This function returns the calculated potential 
    std::vector<double> getPotential();
    /// This function returns the calculated electron concentration
    std::vector<double> getElectronConcentration();
    /// This function returns the calculated hole concentration
    std::vector<double> getHoleConcentration();
    /// This function returns the calculated charge density
    std::vector<double> getChargeDensity();
    ///This function returns the electron-fermi level
    std::vector<double> getEfn();
    ///This function returns the hole-fermi level
    std::vector<double> getEfp();
};





