#include "PSSolver2D.hpp"

PSSolver2D::PSSolver2D(){
    /**
     * Contructor \n 
     * sets dim to 3. dim contains the number of points for one delaunay triangle
     */
    dim = 3;
}

double PSSolver2D::getdistance(long center, long neighbor) {
    /**
     * \n 
     * -center: index of the first point
     * -neighbor: index of the second point
     * \n 
     * \n 
     * The distance between the points with the indices center and neighbor is calculated. The index is referred
     *  to the std::vector points which stores all points
     */
    Vec2d I,J;
    I = points[center];
    J = points[neighbor];
    return std::sqrt((I-J)*(I-J));
}

 double PSSolver2D::getsurface(long center, long neighbor) {
    /**
     * \n 
     * -center: first point's index
     * -neighbor: second point's index
     * \n 
     * \n 
     * First the indices of the triangles which contains both points are searched. The distance between the circumcenter of the triangles is 
     * calculated and returned.
     */
     if(center < neighbor) {
        if(surfaces.find(center) != surfaces.end()) {
            if(surfaces[center].find(neighbor) != surfaces[center].end()) return surfaces[center][neighbor];
        }
    }
    else {
        if(surfaces.find(neighbor) != surfaces.end()) {
            if(surfaces[neighbor].find(center) != surfaces[neighbor].end()) return surfaces[neighbor][center];
        }
    }
     Vec2d I,J;
     long index1 = -1, index2 = -1;
     for(unsigned int i = 0; i < triangles.size(); i++) {
        if(std::find(triangles[i].cbegin(),triangles[i].cend(), center) != triangles[i].cend() && 
        std::find(triangles[i].cbegin(),triangles[i].cend(), neighbor) != triangles[i].cend()) {
            if(index1 == -1) index1 = i;
            else index2 = i;
        }
    }
    if(index2 == -1) {
        J.koord[0] = points[center][0];
        J.koord[1] = points[center][1];
        I.koord[0] = points[neighbor][0];
        I.koord[1] = points[neighbor][1];
        J = (I+J)*0.5;
        I = circumcenter[index1];   
    }
    else {
        I = circumcenter[index1];
        J = circumcenter[index2];
    }
    double result = std::sqrt((I-J)*(I-J));
    if(center < neighbor) {
        surfaces[center][neighbor] = result;
    }
    else {
        surfaces[neighbor][center] = result;
    }
    return result;
 }

 double PSSolver2D::getvolume(long center) {
    /**
     * \n
     * center: the index of the point of which the control volume needs to be calculated
     * \n 
     * \n 
     * First all triangles are searched which contains the center point. The points are stored in the right order and with the class 
     * Surface the control volume is calculated.
     */
    Vec2d A,B, AB;
    std::vector<long> neighbors, fig, copy;
    std::vector<struct List2D> reihenfolge;
    struct List2D new_element;
    double volume = 0;
    std::vector<double> new_point;
    long neighbor_next = -1, neighbor_prev = -1;
    std::fill_n(std::back_inserter(new_point), 2, 0);
    
    if(std::isnan(volumes[center]) == false)  {
        return volumes[center];
    }
       
    //find all triangles that contain center node
    for (unsigned long j = 0; j < triangles.size(); j++) {
      if (triangles[j][0] == center || triangles[j][1] == center || triangles[j][2] == center) {
        fig.push_back(j);
      }
    }

    //add first element to reihenfolge and find both neighbors
   new_element = circumcenter[fig[0]];
   for(long i = 0; i < 3; i++) {
       if(triangles[fig[0]][i] != center) {
           if(new_element.neighbor_next == -1) {
               new_element.neighbor_next = triangles[fig[0]][i];
               neighbor_next = triangles[fig[0]][i];
           }
           else {
               new_element.neighbor_prev = triangles[fig[0]][i];
               neighbor_prev = triangles[fig[0]][i];
           }
       }
   }
   reihenfolge.push_back(new_element);
   fig.erase(fig.begin());

    for(unsigned i = 0; i < fig.size(); i++) {
        for(long t = 0; t < 3; t++) {
            if(triangles[fig[i]][t] == neighbor_next) {
                new_element.neighbor_prev = neighbor_next;
                for(long z = 0; z < 3; z++) {
                    if(triangles[fig[i]][z] != neighbor_next && triangles[fig[i]][z] != center) {
                        neighbor_next = triangles[fig[i]][z];
                        new_element.neighbor_next = neighbor_next;
                        break;
                    }
                }
                new_element = circumcenter[fig[i]];
                reihenfolge.push_back(new_element);
                fig.erase(fig.begin() + i);
                i = -1;
                break;
            }  
        }
    }
    std::reverse(std::begin(reihenfolge), std::end(reihenfolge));
    for(unsigned int i = 0; i < fig.size(); i++) {
        for(long t = 0; t < 3; t++) {
            if(triangles[fig[i]][t] == neighbor_prev) {
                new_element.neighbor_next = neighbor_prev;
                for(long z = 0; z < 3; z++) {
                    if(triangles[fig[i]][z] != neighbor_prev && triangles[fig[i]][z] != center) {
                        neighbor_prev = triangles[fig[i]][z];
                        new_element.neighbor_prev = neighbor_prev;
                        break;
                    }
                }
                new_element = circumcenter[fig[i]];
                reihenfolge.push_back(new_element);
                fig.erase(fig.begin() + i);
                i = -1;
                break;
            }  
        }
    }
    std::reverse(std::begin(reihenfolge), std::end(reihenfolge));
    
    //first and last point are not the same, add buffer points to get a closed loop
    if(reihenfolge[0].neighbor_prev != reihenfolge[reihenfolge.size()-1].neighbor_next) {        
        A = points[center];
        B = points[reihenfolge[reihenfolge.size()-1].neighbor_next];
        AB = (A + B) * 0.5;
        new_element = AB;
        reihenfolge.push_back(new_element);

        A = points[center];
        B = points[reihenfolge[0].neighbor_prev];
        AB = (A + B) * 0.5;
        new_element = AB;

        std::reverse(std::begin(reihenfolge), std::end(reihenfolge));
        reihenfolge.push_back(new_element);
        std::reverse(std::begin(reihenfolge), std::end(reihenfolge));

        new_point[1] = points[center][1];
        new_point[0] = points[center][0];
        new_element = new_point;
        reihenfolge.push_back(new_element);        
    }

    //calculate control volume
    for(unsigned int i = 0; i < reihenfolge.size(); i++) {
        if(i != reihenfolge.size() - 1) {
            volume += reihenfolge[i].p[0] * reihenfolge[i+1].p[1] - reihenfolge[i].p[1] * reihenfolge[i+1].p[0];
        }
        else {
            volume += reihenfolge[i].p[0] * reihenfolge[0].p[1] - reihenfolge[i].p[1] * reihenfolge[0].p[0];
        }
    }
   
    if(volume < 0) volume *= -1;
    volumes[center] = volume/2;
    return volume/2;         
 }


void PSSolver2D::create_mesh() {
    /**
     * This function calculates the circumcenter of every triangle. The circumcenter is calculated with the
     * class Matrix2d.
     */
    std::vector<double> p; 
    std::array<double, 2> right, ergebnis;
    Vec2d I, L, J;
    Matrix2d mat;
    std::fill_n(std::back_inserter(p), 2, 0);

    for(auto var : triangles) {  // edges 1-2, 2-3
        I = points[var[0]];
        J = points[var[1]];
        L = points[var[2]];
         //edge 1-2
        right[0] = ((I + J)*0.5)*(J-I);
        mat.setkoeff(1,1,(J-I).koord[0]);
        mat.setkoeff(1,2,(J-I).koord[1]);
        //edge 2-3
        right[1] = ((L + J)*0.5)*(J-L);
        mat.setkoeff(2,1,(J-L).koord[0]);
        mat.setkoeff(2,2,(J-L).koord[1]);
        ergebnis = mat.solve(right);
        p[0] = ergebnis[0];
        p[1] = ergebnis[1];
        circumcenter.push_back(p);
    } 

    ut = kb*T/q;
    ni = sqrt(Nv*Nc*exp((-(Ec-Ev)/ut)));
    volumes.resize(points.size(), NAN);
   // std::fill_n(std::back_inserter(volumes), indices.size(), 0);
}


void PSSolver2D::import_rectangular_mesh(std::vector<std::array<long, 4> > mesh, long Nx, long Ny) {
    /**
     * \n 
     * -mesh: stores the points indices of the rectangular mesh
     * -Nx:number of points in x direction
     * -Ny: number of points in y direction
     * \n 
     * \n 
     * This function can be used to transform an rectangular mesh into a triangular mesh. I
     */
    (*this).Nx = Nx;
    (*this).Ny = Ny;
    std::vector<long> buffer;
    buffer.resize(3);
    for(unsigned i = 0; i < mesh.size(); i++) {
        buffer[0] = mesh[i][0];
        buffer[1] = mesh[i][1];
        buffer[2] = mesh[i][2];
        triangles.push_back(buffer);

        buffer[0] = mesh[i][2];
        buffer[1] = mesh[i][3];
        buffer[2] = mesh[i][0];
        triangles.push_back(buffer);
    }
    for(unsigned s = 0; s < mesh.size(); s++)  {

        for(auto iterator = semiconductor_segments.begin(); iterator != semiconductor_segments.end(); iterator++) {
            if(std::find((iterator->second).begin(), (iterator->second).end(), s) != end(iterator->second)) {
                for(auto k = mesh[s].begin(); k != mesh[s].end(); k++) { 
                    if(std::find(indices.begin(), indices.end(), *k) == end(indices)) {
                        indices.push_back(*k);
                    }
                }
            }
        }
        for(auto it = contact_segments.begin(); it != contact_segments.end(); it++)  {
            if(std::find(contact_segments.at(it->first).begin(), contact_segments.at(it->first).end(), s) != end(contact_segments.at(it->first))) {
                for(unsigned k = 0; k < 4; k++) {
                    if(std::find(contact_segments_points.at(it->first).begin(), contact_segments_points.at(it->first).end(), mesh[s][k]) == end(contact_segments_points.at(it->first) )) {
                        contact_segments_points.at(it->first).push_back( mesh[s][k]);
                    }
                }
            }
        }
    }
}

void PSSolver2D::set_number_of_points(long Nx, long Ny) {
    (*this).Nx = Nx;
    (*this).Ny = Ny;
}

void PSSolver2D::test() {
    /**
     * \n 
     * This function is only used to make sure the control volumes are calculated
     * 
     */
    double volume = 0;
    create_mesh();
    std::cout << "Nach create mesh! " << std::endl;
    for(unsigned i = 0; i < indices.size(); i++) volume += getvolume(indices[i]);
    std::cout << "Gesamtvolumen " << volume << std::endl;
    std::cout << "indices size " << indices.size() << std::endl;
}

std::vector<double> PSSolver2D::poisson_rect() {
    /**
     * \n 
     * This funcion can be used to solve the poisson equation on a regular rectangular mesh. For regular rectangular meshes 
     * this function should be preferred.   
     */
    double h1, h2, buffer; // buffer, a,b,c,d;
    Eigen::SparseMatrix<double> mat(Nx*Ny,Nx*Ny);
    Eigen::VectorXd vec(Nx*Ny);
    Eigen::VectorXd result(Nx*Ny);

    for(long i = 0; i < Nx*Ny; i++) {
        h1 = 0;
        h2 = 0;
        //find for every Point all Neighbors:
        if((i + 1) % Nx == 0) {
             h1 = getdistance(i, i-1)*0.5;
        }
        else if((i + 1) % Nx == 1) {
            h1 = getdistance(i+1, i)*0.5;
        } 
        else {
            h1 = getdistance(i+1, i-1)*0.5;
        }
        if(i + Nx > (long) Nx*Ny-1) {
            h2 = getdistance(i, i-Nx)*0.5;
        }
        else if (i - Nx < 0){
            h2 = getdistance(i+Nx, i)*0.5;
        }
        else {  
            h2 = getdistance(i+Nx, i-Nx)*0.5;
        }
        vec[i] = 0;
        buffer = 0;

        
        if(dirichlet.find(i) != dirichlet.end()) {
            vec[i] = dirichlet.at(i);
            mat.insert(i,i) = 1;    
            continue;
        }
        
        if(neumann.find(i) != neumann.end()) {
            if((i + 1) % Nx == 0) {
                vec[i] -= h2*neumann.at(i);
            }
            else if((i + 1) % Nx == 1) {
                vec[i] -= h2*neumann.at(i);
            } 
            
            else if(i + Nx > (long) indices.size()-1) {
                vec[i] -= h1*neumann.at(i);
            }
            else if (i - Nx < 0){
                vec[i] -= h1*neumann.at(i);
            }
            else std::cout << "error" << std::endl;
            
        }

        vec[i] -= h1*h2*charge_density[i];
        
        if((i + 1) % Nx == 0) {
            mat.insert(i, i-1) = h2/getdistance(i-1,i);
            buffer -= h2/getdistance(i-1,i);
        }
        else if((i + 1) % Nx == 1) {
            mat.insert(i, i+1) = h2/getdistance(i+1,i);
            buffer -= h2/getdistance(i+1,i);
        } 
        else {
            mat.insert(i, i-1) = h2/getdistance(i-1,i);
            mat.insert(i, i+1) = h2/getdistance(i+1,i);
            buffer -= (h2/getdistance(i-1,i) + h2/getdistance(i+1,i));
        }
        
        if(i + Nx > Nx*Ny -1) {
            mat.insert(i, i-Nx) = h1/getdistance(i-Nx,i);
            buffer -= h1/getdistance(i-Nx,i);
        }
        else if (i - Nx < 0){
            mat.insert(i, i+Nx) = h1/getdistance(i+Nx,i);
            buffer -= h1/getdistance(i+Nx,i);
        }
        else {  
            mat.insert(i, i+Nx) = h1/getdistance(i+Nx,i);
            mat.insert(i, i-Nx) = h1/getdistance(i-Nx,i);
            buffer -= (h1/getdistance(i+Nx,i) + h1/getdistance(i-Nx,i));
        }
                
        mat.insert(i,i) = buffer;
    }

    //std::cout << mat << std::endl << std::endl;
    //std::cout << vec << std::endl << std::endl;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > lscg;
    lscg.compute(mat);
    result = lscg.solve(vec);

    std::vector<double> return_value;
    for(long i = 0; i < result.size(); i++) {
        return_value.push_back(result[i]);
    }
    return return_value;
    
}

void PSSolver2D::scalepoints(std::vector<std::vector<double>> &grid, double factor) {
    for(auto it = grid.begin(); it !=grid.end(); it++) {
        for(auto it2 = it->begin(); it2 != it->end(); it2++) {
            *it2 = *it2*factor;
        }
    }
}


