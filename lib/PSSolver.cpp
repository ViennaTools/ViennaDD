#include "PSSolver.hpp"

void PSSolver::solveDDClassic() {
    /// @brief solves the DD-model. The conservation laws are solved directly for the electron and holeconcentration. Converges faster 
    ///than solveDDQuasiFermi but is also susceptible to unstable behavior. The Gummel method is used to solve the DD-model  \n 
    double max_voltage = 0;
    long n;
    calculateBuiltIn();

    for(auto it = dirichlet.begin(); it != dirichlet.end(); it++) {
        if(std::abs(it->second) > max_voltage) max_voltage = it->second;
    }

    if(max_voltage < std::abs(1e-12)) return; 
    n = static_cast<long>(max_voltage/(ut*steps))+1;
    for(long k = 1; k <= n; k++) { 
        std::cout << "applied voltage " << (max_voltage*k)/n << std::endl;
        if(k == 1) {
            for(auto it = dirichlet.begin(); it != dirichlet.end(); it++) {
                it->second = it->second / n;
            }
        }
        else {
            for(auto it = dirichlet.begin(); it != dirichlet.end(); it++) {
                it->second = it->second *((double)k/(k-1));
            }
        }
        for(long i = 0; i <  max_gummel_iterations; i++) {
            pot = fixpoint(pot);
            elec = electronconcentration(pot);
            holes = holeconcentration(pot);
        } 
    } 
}

void PSSolver::calculateBuiltIn() {
    /// @brief solves the DD-model for thermodynamic equillibrium
    ut = kb*T/q;
    Ei =  ((Ec+Ev)/2)+(ut/2)*log(Nv/Nc);
    ni = sqrt(Nc*Nv*exp(((Ev-Ec)/ut)));
    n1 = ni*exp((Et-Ei)/ut);
    p1 =  ni*exp((Ei-Et)/ut);

    holes.resize(indices.size());
    elec.resize(indices.size());
    pot.resize(indices.size());
    Efn.resize(indices.size());
    Efp.resize(indices.size());
    for(unsigned i = 0; i < indices.size();i++) {
        pot[i] =  (Ec+Ev)/(2*ut) +asinh((dopand[i] - aceptor[i])/(2*ni)) - 0.5*log(Nc/Nv);
        Efn[i] = 0;
        Efp[i] = 0;
    }
    
    intrinsic_debye = sqrt((ut*epsilon_r*8.854*1e-12)/(q*ni*1e6));
    //scale points
    scalepoints(points, 1/intrinsic_debye);
    create_mesh();
    for(long i = 0; i < max_newton_iterations; i++) {
        pot = scaled_nonlinear_poisson(pot);
    }
    for(unsigned i = 0; i < indices.size(); i++) {
        elec[i] = elec_fermi(pot[i]*ut);
        holes[i] = holes_fermi(pot[i]*ut);
        pot[i] = pot[i]*ut;
    }
    scalepoints(points, intrinsic_debye);
}


Eigen::VectorXd PSSolver::scaled_nonlinear_poisson(Eigen::VectorXd potold) {
     /**
     * \n 
     * -oldpot: This is previously claculated potential \n
     * 
     * This function solves a modified poisson equation which is optained by fixedpoint problem. This function is used to caclulate 
     * the new potential with the gummel method.
     */
    Eigen::SparseMatrix<double> mat(indices.size(),indices.size());
    Eigen::VectorXd vec(indices.size());
    Eigen::VectorXd ergebnis(indices.size());
    double buffer;
    std::vector<long> neighbors{};
    long index;
    for(unsigned long i = 0; i < indices.size(); i++) {
         //find for every Point all Neighbors:
        neighbors = find_neighbors(indices[i]);
        vec[i] = 0;
        buffer = 0;
        if(dirichlet.find(indices[i]) != dirichlet.end()) {
            mat.insert(i,i) = 1;
            vec[i] = potold[i];
            continue;
        }
        
        buffer = getvolume(indices[i])*(- exp(-Ei/ut)* exp(potold[i]) - exp(Ei/ut)*exp(-potold[i]));
        vec[i] = (-potold[i]*(exp(-Ei/ut)*exp(potold[i]) + exp(Ei/ut)* exp(-potold[i])) + 
                exp(-Ei/ut)*exp(potold[i]) - exp(Ei/ut)*exp(-potold[i]) -(dopand[i]-aceptor[i])/ni)*getvolume(indices[i]);
            
        for(unsigned long j = 0; j < neighbors.size(); j++) {
            //initialize matrix
            index = std::distance(indices.begin(), find(indices.begin(), indices.end(), neighbors[j]));             
            mat.insert(i,index) = (getsurface(indices[i], neighbors[j])/getdistance(indices[i], neighbors[j]));
            buffer -= (getsurface(indices[i], neighbors[j])/getdistance(indices[i], neighbors[j]));
        }
        
        mat.insert(i,i) = buffer;
        neighbors.clear();
    }

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > lscg;
    lscg.compute(mat);
    ergebnis = lscg.solve(vec);
    return ergebnis;
}


Eigen::VectorXd PSSolver::fixpoint(Eigen::VectorXd oldpot) {
     /**
     * \n 
     * -oldpot: This is previously claculated potential \n
     * 
     * This function solves a modified poisson equation which is optained by fixedpoint problem. This function is used to caclulate 
     * the new potential with the gummel method.
     */

    Eigen::SparseMatrix<double> mat(indices.size(),indices.size());
    Eigen::VectorXd vec(indices.size());
    Eigen::VectorXd ergebnis(indices.size());
    double buffer;
    double  faktor;
    faktor = (1.602176634*1e-1)/(8.854*epsilon_r);
    std::vector<long> neighbors{};
    long index;
    for(unsigned long i = 0; i < indices.size(); i++) {
         //find for every Point all Neighbors:
        neighbors = find_neighbors(indices[i]);
        vec[i] = 0;
        buffer = 0;
        if(dirichlet.find(indices[i]) != dirichlet.end()) {
            vec[i] = (Ec+Ev)/(2) +ut*asinh((dopand[i] - aceptor[i])/(2*ni)) - 0.5*ut*log(Nc/Nv) + dirichlet.at(indices[i]);
            mat.insert(i,i) = 1;
            //buffer = 1;
            neighbors.clear();
            continue;
        }
        buffer -= 1*getvolume(indices[i])*(faktor/ut)*(elec[i]+holes[i]); //
        vec[i] += getvolume(indices[i])*(((elec[i]-dopand[i]-holes[i]+aceptor[i])*1e-1*1.602)/(8.854*epsilon_r)-(faktor/ut)*(elec[i]+holes[i])*oldpot[i]);
            
        for(unsigned long j = 0; j < neighbors.size(); j++) {
            //initialize matrix
            index = std::distance(indices.begin(), find(indices.begin(), indices.end(), neighbors[j]));          
            mat.insert(i,index) = (getsurface(indices[i], neighbors[j])/getdistance(indices[i], neighbors[j]));
            buffer -= (getsurface(indices[i], neighbors[j])/getdistance(indices[i], neighbors[j]));
        }
        mat.insert(i,i) = buffer;
        neighbors.clear();
    }

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > lscg;
    lscg.setMaxIterations(100000);
    lscg.setTolerance(1e-30);
    lscg.compute(mat);
    ergebnis = lscg.solve(vec);
    return ergebnis;
}




std::vector<double> PSSolver::solvePoisson() {
    /**
     * \n 
     * \n 
     * This function solves the classical poisson equation for the stored charge_density.
     */

    create_mesh();
    Eigen::SparseMatrix<double> mat(indices.size(),indices.size());
    Eigen::VectorXd vec(indices.size());
    Eigen::VectorXd ergebnis(indices.size());
    double buffer;
    std::vector<long> neighbors{};
    long index;
    sort(indices.begin(), indices.end());
    for(unsigned long i = 0; i < indices.size(); i++) {
         //find for every Point all Neighbors:
        neighbors = find_neighbors(indices[i]);
        vec[i] = 0;
        buffer = 0;
        if(dirichlet.find(indices[i]) != dirichlet.end()) {
            vec[i] = dirichlet.at(indices[i]);
            mat.insert(i,i) = 1;
            neighbors.clear();
            continue;
        }
        
        if(neumann.find(indices[i]) != neumann.end()) {
            long counter;
            for(auto it = neighbors.begin(); it != neighbors.end(); it++) {
                counter = 0;
                for(auto it2 = triangles.begin(); it2 != triangles.end(); it2++) {
                    if(std::find(it2->begin(), it2->end(), *it) != it2->end() &&
                        std::find(it2->begin(), it2->end(), indices[i]) != it2->end()) counter++;
                }
                if(counter == 1) {
                    vec[i] -= getdistance(indices[i], *it)*0.5*neumann.at(indices[i]);      //*1e+12/(8.854*epsilon_r);
                }
            }
        }
        vec[i] -= getvolume(indices[i])*charge_density[i];
        for(unsigned long j = 0; j < neighbors.size(); j++) {
            //initialize matrix
            index = std::distance(indices.begin(), std::find(indices.begin(), indices.end(), neighbors[j]));          
            mat.insert(i,index) = (getsurface(indices[i], neighbors[j])/getdistance(indices[i], neighbors[j]));
            buffer -= (getsurface(indices[i], neighbors[j])/getdistance(indices[i], neighbors[j]));
        }
        mat.insert(i,i) = buffer;
        neighbors.clear();
    }
    //loesen
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > lscg;
    lscg.compute(mat);
    ergebnis = lscg.solve(vec);
    std::vector<double> res;
    for(int i = 0; i < ergebnis.size(); i++) res.push_back(ergebnis[i]);
    return res;
}

Eigen::VectorXd PSSolver::electronconcentration(Eigen::VectorXd potential) {
    /**
     * \n
     * potential: is the current potential for which the electronconcentration is calulated
     * \n \n 
     * This function calculates the electronconcentration for the current potential with the Scharfetter-Gummel-scheme. The unit of the 
     * returned electron concentration is 1/cm³.
     */
    Eigen::VectorXd vec(indices.size());
    Eigen::SparseMatrix<double> mat(indices.size(), indices.size());
    std::vector<long> neighbors{};
    double buffer;
    long index;
    
    //Scharfetter-gummel-scheme
    for(unsigned long i = 0; i < indices.size(); i++) {
        //find for every Point all Neighbors:
        neighbors = find_neighbors(indices[i]);
        
        if(std::find(contact_segments_points.at("Anode").begin(), contact_segments_points.at("Anode").end(), indices[i]) != end(contact_segments_points.at("Anode"))) {
            vec[i] = (0.5*(dopand[i]-aceptor[i])+sqrt(0.25*pow(dopand[i]-aceptor[i],2)+pow(ni,2)))*1e6;
            mat.insert(i,i) =  1;
        }
        else if(std::find(contact_segments_points.at("Kathode").begin(), contact_segments_points.at("Kathode").end(), indices[i]) != end(contact_segments_points.at("Kathode"))) {
            vec[i] = (0.5*(dopand[i]-aceptor[i])+sqrt(0.25*pow(dopand[i]-aceptor[i],2)+pow(ni,2)))*1e6;
            mat.insert(i,i) =  1;
        }
        else {
            buffer = 0;
            vec[i] =0;
            for(unsigned long j = 0; j < neighbors.size(); j++) {
                //initialize matrix
                index = std::distance(indices.begin(), std::find(indices.begin(), indices.end(), neighbors[j]));          
                mat.insert(i,index) = (1/getdistance(indices[i], neighbors[j]))
                        *bernoulifunction(potential[i],potential[index])*getsurface(indices[i], neighbors[j]);
                buffer -=(1/getdistance(indices[i], neighbors[j]))
                        *bernoulifunction(potential[index],potential[i])*getsurface(indices[i], neighbors[j]);
                }
            mat.insert(i,i) = buffer;
        }
        neighbors.clear();
    }

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> lscg;
    lscg.setMaxIterations(100000);
    lscg.setTolerance(1e-30);
    lscg.compute(mat);
    Eigen::VectorXd ergebnis = lscg.solve(vec);
    for(long i = 0; i < ergebnis.size(); i++) ergebnis[i] = ergebnis[i] /1e6;
    return ergebnis;

}

Eigen::VectorXd PSSolver::holeconcentration(Eigen::VectorXd potential) {
    /**
     * \n
     * potential: is the current potential for which the electronconcentration is calulated
     * \n \n 
     * This function calculates the hole concentration for the current potential with the Scharfetter-Gummel-scheme. The unit of the 
     * returned hole concentration is 1/cm³.
     */
    Eigen::VectorXd  vec(indices.size());
    Eigen::SparseMatrix<double> mat(indices.size(), indices.size());
    std::vector<long> neighbors{};
    double buffer;
    long index;
   
    for(unsigned long i = 0; i < indices.size(); i++) {
         //find for every Point all Neighbors:
        neighbors = find_neighbors(indices[i]);
        
        if(std::find(contact_segments_points.at("Anode").begin(), contact_segments_points.at("Anode").end(), indices[i]) != end(contact_segments_points.at("Anode"))) {
            vec[i] = (-0.5*(dopand[i]-aceptor[i])+sqrt(0.25*pow(dopand[i]-aceptor[i],2)+pow(ni,2)))*1e6;
            mat.insert(i,i) =  1;
        }
        else if(std::find(contact_segments_points.at("Kathode").begin(), contact_segments_points.at("Kathode").end(), indices[i]) != end(contact_segments_points.at("Kathode"))) {
            vec[i] = (-0.5*(dopand[i]-aceptor[i])+sqrt(0.25*pow(dopand[i]-aceptor[i],2)+pow(ni,2)))*1e6;
            mat.insert(i,i) =  1;
        }
        else {
            buffer = 0;
            vec[i] = 0;
            for(unsigned long j = 0; j < neighbors.size(); j++) {
                //initialize matrix
                index = std::distance(indices.begin(), std::find(indices.begin(), indices.end(), neighbors[j]));          
                mat.insert(i,index) = (-1/getdistance(indices[i], neighbors[j]))
                        *bernoulifunction(potential[index],potential[i])*getsurface(indices[i], neighbors[j]);
                buffer +=(1/getdistance(indices[i], neighbors[j]))
                        *bernoulifunction(potential[i],potential[index])*getsurface(indices[i], neighbors[j]);
                }
            mat.insert(i,i) = buffer;
        }
        neighbors.clear();
    }

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> lscg;
    lscg.setMaxIterations(10000);
    lscg.setTolerance(1e-30);
    lscg.compute(mat);
    Eigen::VectorXd ergebnis = lscg.solve(vec);
    for(long i = 0; i < ergebnis.size(); i++) ergebnis[i] = ergebnis[i] /1e6;
    return ergebnis;
}

Eigen::VectorXd PSSolver::newtoniteration(Eigen::VectorXd potold) {
    /**
     * \n
     * -potold: is the previously calculated potential
     * \n 
     * \n 
     * This function solves the nonlinear poisson equation in equilibrium. This function perfoms one newtoniteration to 
     * calculate the delta-potential. This delta-potential is then added to potold and the the result is returned.
     */
    Eigen::SparseMatrix<double> mat(indices.size(),indices.size());
    Eigen::VectorXd vec(indices.size());
    Eigen::VectorXd ergebnis(indices.size());
    double buffer;
    double  faktor;
    faktor = (1.602176634*1e-1)/(8.854*epsilon_r);
    std::vector<long> neighbors{};
    long index;

    for(unsigned long i = 0; i < indices.size(); i++) {
         //find for every Point all Neighbors:
        for(unsigned long j = 0; j < triangles.size(); j++) {
            if(std::find(triangles[j].begin(), triangles[j].end(), indices[i]) != triangles[j].end()) {
               for(long k = 0; k < dim; k++) {
                   if(triangles[j][k] != indices[i] && 
                    std::find(neighbors.begin(), neighbors.end(), triangles[j][k]) == neighbors.end()) {
                        neighbors.push_back(triangles[j][k]);
                    }
               }
            }
        }
        if(std::find(contact_segments_points.at("Anode").begin(), contact_segments_points.at("Anode").end(), indices[i]) != end(contact_segments_points.at("Anode"))) {
            vec[i] = -1*(faktor*(-elec_fermi(potold[i])+holes_fermi(potold[i])+dopand[i]-aceptor[i]));
            mat.insert(i,i) =  - (faktor)*(elec_fermi(potold[i])+holes_fermi(potold[i]))/ut;
        }
        else if(std::find(contact_segments_points.at("Kathode").begin(), contact_segments_points.at("Kathode").end(), indices[i]) != end(contact_segments_points.at("Kathode"))) {
            vec[i] = -1*(faktor*(-elec_fermi(potold[i])+holes_fermi(potold[i])+dopand[i]-aceptor[i]));
            mat.insert(i,i) =  - (faktor)*(elec_fermi(potold[i])+holes_fermi(potold[i]))/ut;
        }
        else {

            buffer = 0;
            buffer = ((getvolume(indices[i])*faktor)/ut)*(-elec_fermi(potold[i])-holes_fermi(potold[i]));
            vec[i] = getvolume(indices[i])*((elec_fermi(potold[i])-holes_fermi(potold[i])-dopand[i]+aceptor[i])*1e-1*1.602)/(8.854*epsilon_r);
            for(unsigned long j = 0; j < neighbors.size(); j++) {
                //initialize matrix
                index = std::distance(indices.begin(), std::find(indices.begin(), indices.end(), neighbors[j]));          
                vec[i] -= (getsurface(indices[i], neighbors[j])/getdistance(indices[i], neighbors[j]))*(potold[index]-potold[i]);
                buffer -= getsurface(indices[i], neighbors[j])/getdistance(indices[i], neighbors[j]);
                
                mat.insert(i,index) = getsurface(indices[i], neighbors[j])/getdistance(indices[i], neighbors[j]);
            }
            mat.insert(i,i) = buffer; 
        }
        neighbors.clear();
    }

    //loesen
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > lscg;
    lscg.setMaxIterations(100000);
    lscg.setTolerance(1e-30);
    lscg.compute(mat);
    ergebnis = lscg.solve(vec);
    for(long i = 0; i < ergebnis.size(); i++) {
        ergebnis[i] = ergebnis[i]+ potold[i];
    }
  
   
    return ergebnis;
}


double PSSolver::error(Eigen::VectorXd potold) {
     /**
     * \n 
     * -potential: thi is the current potential
     * \n 
     * \n 
     * The maximum  error of newton's method is calculated and returned. This error is used as a break criteria 
     */
    double buffer;
    double error = 0;
    std::vector<long> neighbors{};
    long index;
    for(unsigned long i = 0; i < indices.size(); i++) {
         //find for every Point all Neighbors:
        for(unsigned long j = 0; j < triangles.size(); j++) {
            if(std::find(triangles[j].begin(), triangles[j].end(), indices[i]) != triangles[j].end()) {
               for(long k = 0; k < dim; k++) {
                   if(triangles[j][k] != indices[i] && 
                    std::find(neighbors.begin(), neighbors.end(), triangles[j][k]) == neighbors.end()) {
                        neighbors.push_back(triangles[j][k]);
                    }
               }
            }
        }
        
        if(std::find(contact_segments_points.at("Anode").begin(), contact_segments_points.at("Anode").end(), indices[i]) != end(contact_segments_points.at("Anode"))) {
        }
        else if(std::find(contact_segments_points.at("Kathode").begin(), contact_segments_points.at("Kathode").end(), indices[i]) != end(contact_segments_points.at("Kathode"))) {
        }
        else {
            buffer = 0;
            buffer = -getvolume(indices[i])*((elec_fermi(potold[i])-holes_fermi(potold[i])-dopand[i])*1e-1*1.602)/(8.854*epsilon_r);
            for(unsigned long j = 0; j < neighbors.size(); j++) {
                //initialize matrix
                index = std::distance(indices.begin(), std::find(indices.begin(), indices.end(), neighbors[j]));          
                buffer += (getsurface(indices[i], neighbors[j])/getdistance(indices[i], neighbors[j]))*(potold[index]-potold[i]);
            } 
        }
        error += std::abs(buffer);
        neighbors.clear();
    }
    return error;
}



void PSSolver::setCellGrid(std::vector<std::vector<long> > data) {
    /**
     * \n 
     * -data: containes the indices of the points which are connected.
     * \n 
     * \n 
     * This function has to be called after addSemiconductorSegment(), addContactSegment() and setSimulationGrid().
     * 
     */
    for(unsigned long s = 0;s < data.size(); s++) {
        for(auto iterator = semiconductor_segments.begin(); iterator != semiconductor_segments.end(); iterator++) {
            if(std::find((iterator->second).begin(), (iterator->second).end(), s) != end(iterator->second)) {
                triangles.push_back(data[s]);
                for(auto k = data[s].begin(); k != data[s].end(); k++) { 
                    if(std::find(indices.begin(), indices.end(), *k) == end(indices)) {
                        indices.push_back(*k);
                    }
                }
            }
        }
        auto it =  contact_segments.begin();
        while(it != contact_segments.end()) {
            if(std::find(contact_segments.at(it->first).begin(), contact_segments.at(it->first).end(), s) != end(contact_segments.at(it->first))) {
                for(unsigned k = 0; k < dim; k++) {
                    if(std::find(contact_segments_points.at(it->first).begin(), contact_segments_points.at(it->first).end(), data[s][k]) == end(contact_segments_points.at(it->first) )) {
                        contact_segments_points.at(it->first).push_back( data[s][k]);
                    }
                }
            }
            it++;
        }
    }
}

void PSSolver::setSimulationGrid(std::vector<std::vector<double> > grid) {
    points = grid;
}

 void PSSolver::setDonorConcentration(std::vector<double>  concentration) {
    /**
     * \n 
     * -concentration: donor concentration in 1/cm³. is defined on each point of the simulation grid
     * \n 
     * \n 
     * This function has to be called after setCellGrid().
     */
    if(dopand.size() != indices.size()) {
        dopand.resize(indices.size());
        aceptor.resize(indices.size());
    }
    for(unsigned long i = 0; i < indices.size(); i++) {
        if(isnan(concentration[indices[i]])== 0) {
            dopand[i] = concentration[indices[i]];
        }
    }
}

void PSSolver::setAcceptorConcentration(std::vector<double>  concentration) {
    /**
     * \n 
     * -concentration: aceptor concentration in 1/cm³. is defined on each point of the simulation grid
     * \n 
     * \n 
     * This function has to be called after setCellGrid().
     */
    if(dopand.size() != indices.size()) {
        dopand.resize(indices.size());
        aceptor.resize(indices.size());
    }
    for(unsigned long i = 0; i < indices.size(); i++) {
        if(isnan(concentration[indices[i]])== 0) {
            aceptor[i] = concentration[indices[i]];
        }
    }
}

double PSSolver::bernoulifunction(double Vi, double Vj) {
    double parameter;
    parameter = (Vj-Vi)/ut;
    //without exp   
    if(std::abs(parameter)<=1e-24) return 1 - parameter/2;
    else return parameter/(exp(parameter)-1);
}  

void PSSolver::setNumberGummelIterations(long num) {
    max_gummel_iterations = num;
}

void PSSolver::setNumberNewtonIterations(long num) {
    max_newton_iterations = num;
}

void PSSolver::addContactSegment(std::vector<long> segment, std::string name) {
    /**
     * \n 
     * -name: the name of the contact segment which is added
     * -segment: contains the indices of the triangles which are added to the contact segment
     * \n 
     * \n 
     * The indices of the triangles which are added as a contact segment are stored in a std::map. This funtion must be called before setCellGrid().
     */

    contact_segments.insert(std::make_pair(name, segment));
    contact_segments_points.insert(std::make_pair(name,0));
}

void PSSolver::addContactPoints(std::vector<long> points, std::string name) {
    /**
     * \n 
     * -name: the name of the contact segment which is added
     * -points: The indices of the points which are added to the contact segment name.
     * \n 
     * \n 
     * The indices of the points which are added to a contact segment are stored in a std::map. This funtion can be called after setCellGrid().
     */
    contact_segments_points.insert(std::make_pair(name, points));
}

void PSSolver::addSemiconductorSegment(std::vector<long> data, std::string name) {
    /**
     * \n 
     * -name: name of the semiconductor segment which is added
     * -data: The indices of the triangles which are added to the semiconductor segment name.
     * \n 
     * \n 
     * The indices of the triangles which are added to a semiconductor segment are stored in a std::map. This funtion has to be called before setCellGrid().
     */
    semiconductor_segments.insert(std::make_pair(name, data));
}

void PSSolver::setEpsilonGummel(double eps) {
    epsilongummel = eps;
}

void PSSolver::setPermitivity(double perm) {
    epsilon_r = perm;
}

void PSSolver::setTemperatur(double temp) {
    T = temp;
    ut = (T*1.3806*1e-4)/1.602;
}

double PSSolver::holes_fermi(double voltage) {
    /**
     * -voltage: current potential
     * \n
     * This function returns the hole concentration for the current voltage. The boltzman statistics is used
     */
    return ni*exp(Ei/ut)*exp(-voltage/ut);
    //return Nv*exp((Ev-voltage)/ut);
}

double PSSolver::elec_fermi(double voltage) {
    /**
     * -voltage: current potential
     * \n
     * This function returns the electron concentration for the current voltage. The boltzman statistics is used
     */
    return ni*exp(-Ei/ut)* exp(voltage/ut);
    //return Nc*exp((voltage-Ec)/ut);
}

void PSSolver::setEpsilonNewton(double eps) {
    epsilonnewton = eps;
}

void PSSolver::setElecDensityofStates(double value) {
    Nc = value;
}
void PSSolver::setHolesDensityofStates(double value) {
    Nv = value;
}

void PSSolver::applyVoltageSegment(double voltage, std::string name) {
    /**
     * \n 
     * -voltage: This voltage is applied to the contact segment name
     * -name: name of the contact segment. This segment must exist
     * \n
     * \n
     * If the segment is created with addContactSegment() the function has to be called after setCellGrid(). If the segment is created with 
     * addContactPoints() it can be called before setCellGrid().
     */
    if(contact_segments_points.find(name) != contact_segments_points.end()) {
        auto buffer = contact_segments_points.at(name);
        for(unsigned i = 0; i < buffer.size(); i++) {
            if(dirichlet.find(buffer[i]) == dirichlet.end()) {
                dirichlet.insert(std::make_pair(buffer[i], voltage));
            }
            else {
                dirichlet[buffer[i]] = voltage;
            }
        }
    }
    else {
        std::cout << "Segment " << name << " not found" << std::endl;
    }
    
}


void PSSolver::deleteVoltageSegment(std::string name) {
    /***
     * -name: The applied voltage to this segment is deleted
    */
    if(contact_segments_points.find(name) != contact_segments_points.end()) {
        auto buffer = contact_segments_points.at(name);
        for(unsigned i = 0; i < buffer.size(); i++) {
            if(dirichlet.find(buffer[i]) != dirichlet.end()) {
                dirichlet.erase(buffer[i]);
            }
        }
    }
    else {
        std::cout << "Segment " << name << " not found" << std::endl;
    }
}

void PSSolver::applyVoltagePoint(double voltage, long index){
    if(dirichlet.find(index) == dirichlet.end()) {
        dirichlet.insert(std::make_pair(index, voltage));
    }
    else {
        dirichlet[index] = voltage;
    }
}

void PSSolver::deleteVoltagePoint(long index) {
    if(dirichlet.find(index) != dirichlet.end()) {
        dirichlet.erase(index);
    }
}


void PSSolver::addNeumannSegment(double value, std::string name) {
    /**
     * \n 
     * -value: is applied as a neumann boundary condition to the segment name
     * -name: name of the contact segment. This segment must exist
     * \n
     * \n
     * If the segment is created with addContactSegment() the function must be called after setCellGrid(). If the segment is created with 
     * addContactPoints() it can be called before setCellGrid().
     */
    if(contact_segments_points.find(name) != contact_segments_points.end()) {
        auto buffer = contact_segments_points.at(name);
        for(unsigned i = 0; i < buffer.size(); i++) {
            if(neumann.find(buffer[i]) == neumann.end()) {
                neumann.insert(std::make_pair(buffer[i], value));
            }
            else {
                neumann[buffer[i]] = value;
            }
        }
    }
    else {
        std::cout << "Segment " << name << " not found" << std::endl;
    }
}

void PSSolver::deleteNeumannSegment(std::string name) {
    if(contact_segments_points.find(name) != contact_segments_points.end()) {
        auto buffer = contact_segments_points.at(name);
        for(unsigned i = 0; i < buffer.size(); i++) {
            if(neumann.find(buffer[i]) != neumann.end()) {
                neumann.erase(buffer[i]);
            }
        }
    }
    else {
        std::cout << "Segment " << name << " not found" << std::endl;
    }
}

void PSSolver::addNeumannPoint(double value, long index) {
     if(neumann.find(index) == neumann.end()) {
        neumann.insert(std::make_pair(index, value));
    }
    else {
        neumann[index] = value;
    }
}

void PSSolver::deleteNeumannPoint(long index) {
    if(neumann.find(index) != neumann.end()) {
        neumann.erase(index);
    }
}

void PSSolver::setChargeDensity(std::vector<double> charge) {
    charge_density = charge;
}

void PSSolver::setConductionBandEnergy(double _Ec) {
    Ec = _Ec;
}
void PSSolver::setValenceBandEnergy(double _Ev) {
    Ev = _Ev;
}

void PSSolver::create_mesh(){ 
}
double PSSolver::getsurface(long center, long neighbor){ return -1*center*neighbor;}
double PSSolver::getvolume(long center){ return -1*center;}
double PSSolver::getdistance(long center, long neighbor){return -1*center*neighbor;}


std::vector<long> PSSolver::find_neighbors(long center) {
    /***
     * returns the indices of all neighbors to the center node.
    */
    std::vector<long> neighbors;
    for(unsigned long j = 0; j < triangles.size(); j++) {
        if(std::find(triangles[j].begin(), triangles[j].end(), center) != triangles[j].end()) {
            for(long k = 0; k < dim; k++) {
                if(triangles[j][k] != center && 
                std::find(neighbors.begin(), neighbors.end(), triangles[j][k]) == neighbors.end()) {
                    neighbors.push_back(triangles[j][k]);
                }
            }
        }
    }
    return neighbors;
}


std::vector<double> PSSolver::getPotential() {
    std::vector<double> result;
    result.resize(points.size(), NAN);
    for(unsigned i = 0; i < indices.size(); i++) {
        result[indices[i]] = pot[i]; 
    }
    return result;
}

std::vector<double> PSSolver::getElectronConcentration() {
    std::vector<double> result;
    result.resize(points.size(), NAN);
    for(unsigned i = 0; i < indices.size(); i++) {
        result[indices[i]] = elec[i]; 
    }
    return result;
}

std::vector<double> PSSolver::getHoleConcentration() {
    std::vector<double> result;
    result.resize(points.size(), NAN);
    for(unsigned i = 0; i < indices.size(); i++) {
        result[indices[i]] = holes[i]; 
    }
    return result;
}

std::vector<double> PSSolver::getChargeDensity() {
    std::vector<double> result;
    result.resize(points.size(), NAN);
    for(unsigned i = 0; i < indices.size(); i++) {
        result[indices[i]] = (-elec[i]+holes[i]-aceptor[i]+dopand[i])*q*1e6;
    }
    return result;
}

void PSSolver::scalepoints(std::vector<std::vector<double>> &grid, double factor) {
    for(auto it = grid.begin(); it !=grid.end(); it++) {
        for(auto it2 = it->begin(); it2 != it->end(); it2++) {
            *it2 = *it2*factor;
        }
    }
}


void PSSolver::setSHR(bool status) {
    SHR = status;
}
void PSSolver::setAuger(bool status) {
    Auger = status;
}
void PSSolver::setSpond(bool status) {
    Spond = status;
}

void PSSolver::setElectronMobility(double mob) {
    elec_mobility = mob;
}
void PSSolver::setHoleMobility(double mob) {
    hole_mobility = mob;
}
void PSSolver::setSpontanRecomb(double _rspond) {
    rspond = _rspond;
}
void PSSolver::setCn(double _Cn) {
    Cn = _Cn;
}
void PSSolver::setCp(double _Cp) {
    Cp = _Cp;
}

void PSSolver::setTrapEnergy(double energy) {
    Et = energy;
}

void PSSolver::setElecLifetime(double time) {
    elec_lifetime = time;
}

void PSSolver::setHoleLifetime(double time) {
    hole_lifetime = time;
}


void PSSolver::electronconcentrationfermi() {
    /***
     * Solves the electron concentration by using the quasi-Fermi levels. To solve the system of nonlinear equations Newton's method is used. 
     * The quasi-fermi levels are scaled by the thermal voltage. 
    */
    Eigen::SparseMatrix<double> mat(indices.size(),indices.size());
    Eigen::VectorXd vec(indices.size());
    Eigen::VectorXd ergebnis(indices.size());
    double buffer;
    double volume, surface, distance;
    double n_fac, p_fac, n1_fac, p1_fac;
    std::vector<long> neighbors{};
    long index;
    for(unsigned long i = 0; i < indices.size(); i++) {
         //find for every Point all Neighbors:
        neighbors = find_neighbors(indices[i]);
        vec[i] = 0;
        buffer = 0;
        if(dirichlet.find(indices[i]) != dirichlet.end()) {
            mat.insert(i,i) = 1;
            vec[i] = -Efn[i]+dirichlet.at(indices[i])/ut; 
            continue;
        }
        volume = getvolume(indices[i]);
        n_fac = exp(-Efn[i])*exp((pot[i]-Ei)/ut);
        p_fac = exp(Efp[i])*exp((Ei-pot[i])/ut);
        n1_fac = exp((Ei-Et)/ut);
        p1_fac = exp((Et-Ei)/ut);
        buffer = 0;
        vec[i] = 0;

        if(Spond == true) {
            vec[i] += volume*ni*rspond*(exp(Efp[i]-Efn[i])-1);     
            buffer += volume*ni*rspond*(exp(Efp[i]-Efn[i]))/ut;
        }
        else if(Auger == true) {
            vec[i] += volume*ni*ni*(Cn*exp((pot[i]-Ei)/ut)*exp(-Efn[i]) + Cp*exp((Ei-pot[i])/ut)*exp(Efp[i]))*(exp(Efp[i]-Efn[i])-1);
        
            buffer += volume*ni*ni*(Cn*exp((pot[i]-Ei)/ut)*exp(-Efn[i]) + Cp*exp((Ei-pot[i])/ut)*exp(Efp[i]))*(exp(Efp[i]-Efn[i]))/ut +
                            volume*ni*ni*(Cn*exp((pot[i]-Ei)/ut)*exp(-Efn[i]))*(exp(Efp[i]-Efn[i])-1)/ut;        
        }
        else if(SHR == true) {
            vec[i] += volume*(((exp(Efp[i]-Efn[i])-1))/(hole_lifetime*(n_fac+n1_fac)+elec_lifetime*(p_fac+p1_fac)));
              
            buffer += volume* (1/ut)*((exp(Efp[i]-Efn[i])*(hole_lifetime*(n_fac+n1_fac)+elec_lifetime*(p_fac+p1_fac))-(exp(Efp[i]-Efn[i])-1)*hole_lifetime*n_fac)/(pow((elec_lifetime*(p_fac+p1_fac)+hole_lifetime*(n_fac+n1_fac)),2)));
        }
        
        buffer = 0;
        vec[i] = 0;
        for(unsigned long j = 0; j < neighbors.size(); j++) {
            
            index = std::distance(indices.begin(), find(indices.begin(), indices.end(), neighbors[j]));
            surface = getsurface(indices[i], neighbors[j]);
            distance = getdistance(indices[i], neighbors[j]);
            vec[i] -= ((surface*ut*elec_mobility*1e-4)/distance)*(exp(-Efn[index])*exp((pot[index]-Ei)/ut)*bernoulifunction(pot[i],pot[index]) - exp(-Efn[i])*exp((pot[i]-Ei)/ut)*bernoulifunction(pot[index],pot[i])); 
            
            mat.insert(i,index) = -((surface*elec_mobility*1e-4)/distance)*exp(-Efn[index])*exp((pot[index]-Ei)/ut)*bernoulifunction(pot[i],pot[index]);
            buffer += ((surface*elec_mobility*1e-4)/distance)*(exp(-Efn[i])*exp((pot[i]-Ei)/ut)*bernoulifunction(pot[index],pot[i]));    
        }
        
        mat.insert(i,i) = buffer;
        neighbors.clear();
    }

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > lscg;
    lscg.setMaxIterations(indices.size());
    lscg.setTolerance(1e-30);
    lscg.compute(mat);
    ergebnis = lscg.solve(vec);
    for(unsigned i = 0; i < indices.size(); i++) {
        Efn[i] = Efn[i] + ergebnis[i];
        elec[i] = ni*exp((pot[i]-Ei)/ut)*exp(-Efn[i]);
    }
}



void PSSolver::holeconcentrationfermi() {
    /***
     * Solves the hole concentration by using the quasi-Fermi levels. To solve the system of nonlinear equations Newton's method is used. 
     * The quasi-fermi levels are scaled by the thermal voltage. 
    */
    Eigen::SparseMatrix<double> mat(indices.size(),indices.size());
    Eigen::VectorXd vec(indices.size());
    Eigen::VectorXd ergebnis(indices.size());
    double buffer;
    double volume, surface, distance;
    double n_fac, p_fac, n1_fac, p1_fac;
    std::vector<long> neighbors{};
    long index;
    for(unsigned long i = 0; i < indices.size(); i++) {
         //find for every Point all Neighbors:
        neighbors = find_neighbors(indices[i]);
        if(dirichlet.find(indices[i]) != dirichlet.end()) {
            mat.insert(i,i) = 1;
            vec[i] = -Efp[i]+dirichlet.at(indices[i])/ut; 
            continue;
        }
        volume = getvolume(indices[i]);
        n_fac = exp(-Efn[i])*exp((pot[i]-Ei)/ut);
        p_fac = exp(Efp[i])*exp((Ei-pot[i])/ut);
        n1_fac = exp((Ei-Et)/ut);
        p1_fac = exp((Et-Ei)/ut);
        buffer = 0;
        vec[i] = 0;

       if(Spond == true) {
            vec[i] -= volume*ni*rspond*(exp(Efp[i]-Efn[i])-1);
            buffer += volume*ni*rspond*(exp(Efp[i]-Efn[i]))/ut;
              
        }
        else if(Auger == true) {
            vec[i] -= volume*ni*ni*(Cn*exp((pot[i]-Ei)/ut)*exp(-Efn[i]) + Cp*exp((Ei-pot[i])/ut)*exp(Efp[i]))*(exp(Efp[i]-Efn[i])-1);
            buffer += volume*ni*ni*(Cn*exp((pot[i]-Ei)/ut)*exp(-Efn[i]) + Cp*exp((Ei-pot[i])/ut)*exp(Efp[i]))*(exp(Efp[i]-Efn[i]))/ut +
                        volume*ni*ni*(Cp*exp((Ei-pot[i])/ut)*exp(Efp[i]))*(exp(Efp[i]-Efn[i])-1)/ut;
                     
        }
        else if(SHR == true) {
            vec[i] -= volume*(((exp(Efp[i]-Efn[i])-1))/(hole_lifetime*(n_fac+n1_fac)+elec_lifetime*(p_fac+p1_fac)));
                buffer += volume* (1/ut)*((exp(Efp[i]-Efn[i])*(hole_lifetime*(n_fac+n1_fac)+elec_lifetime*(p_fac+p1_fac))-(exp(Efp[i]-Efn[i])-1)*elec_lifetime*p_fac)/
                    (pow((elec_lifetime*(p_fac+p1_fac)+hole_lifetime*(n_fac+n1_fac)),2)));
        }
        
        buffer = 0;
        vec[i] = 0;
        for(unsigned long j = 0; j < neighbors.size(); j++) {
            
            index = std::distance(indices.begin(), find(indices.begin(), indices.end(), neighbors[j]));
            surface = getsurface(indices[i], neighbors[j]);
            distance = getdistance(indices[i], neighbors[j]);
            
            
            vec[i] +=  ((surface*ut*hole_mobility*1e-4)/distance)*(exp(Efp[index])*exp((Ei-pot[index])/ut)*bernoulifunction(pot[index],pot[i]) - exp(Efp[i])*exp((Ei-pot[i])/ut)*bernoulifunction(pot[i],pot[index]));
            
            mat.insert(i,index) = -((surface*hole_mobility*1e-4)/distance)*exp(Efp[index])*exp((Ei-pot[index])/ut)*bernoulifunction(pot[index],pot[i]);
            buffer += ((surface*hole_mobility*1e-4)/distance) * exp(Efp[i])*exp((Ei-pot[i])/ut)*bernoulifunction(pot[i],pot[index]);  
                   
        }
        
        mat.insert(i,i) = buffer;
        neighbors.clear();
    }

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > lscg;
    lscg.setMaxIterations(indices.size());
    lscg.setTolerance(1e-30);
    lscg.compute(mat);
    ergebnis = lscg.solve(vec);
    for(unsigned i = 0; i < indices.size(); i++) {
        Efp[i] = Efp[i] +ergebnis[i];
        holes[i] = ni*exp((Ei-pot[i])/ut)*exp(Efp[i]);
    }
}


void PSSolver::solveDDQuasiFermi() {
    // @brief solves the DD-model. The conservation laws are solved for the quasi-fermi levels of electrons and holes. Converges slower 
    ///than solveDDClassic but has a more stable behaviour. The Gummel method is used to solve the DD-model \n 
    double max_voltage = 0;
    long n;
    calculateBuiltIn();

    for(auto it = dirichlet.begin(); it != dirichlet.end(); it++) {
        if(std::abs(it->second) > max_voltage) max_voltage = it->second;
    }
    if(max_voltage < std::abs(1e-12)) return; 
    n = static_cast<long>(max_voltage/(ut*steps))+1;
    for(long k = 1; k <= n; k++) { 
        std::cout << "applied voltage " << (max_voltage*k)/n << std::endl;
        if(k == 1) {
            for(auto it = dirichlet.begin(); it != dirichlet.end(); it++) {
                it->second = it->second / n;
            }
        }
        else {
            for(auto it = dirichlet.begin(); it != dirichlet.end(); it++) {
                it->second = it->second *((double)k/(k-1));
            }
        }
        for(long i = 0; i <  max_gummel_iterations; i++) {
            pot = fixpoint(pot);
            electronconcentrationfermi();
            holeconcentrationfermi();
        } 
    }  
 }


std::vector<double> PSSolver::getEfn() {
    std::vector<double> result;
    for(unsigned i = 0; i < indices.size(); i++) result.push_back(Efn[i]);
    return result;
}
std::vector<double> PSSolver::getEfp() {
    std::vector<double> result;
    for(unsigned i = 0; i < indices.size(); i++) result.push_back(Efp[i]);
    return result;
}

void PSSolver::setSteps(double _steps) {
    steps = _steps;
}
