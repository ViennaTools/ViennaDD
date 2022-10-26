#include "PSSolver3D.hpp"

PSSolver3D::PSSolver3D() {
    /**
     * Contructor \n 
     * sets dim to 3. dim contains the number of points for one delaunay triangle
     */
    dim = 4;
}

void PSSolver3D::create_mesh() {
        /**
     * This function calculates the circumcenter of every triangle. The circumcenter is calculated with the
     * class Matrix3d.
     */
    std::vector<double> p;
    std::array<double, 3> right, ergebnis;
    Vec3d I, L, J, K;
    Matrix3d mat;
    std::fill_n(std::back_inserter(p), 3, 0);

    //calculate circumecenter for all triangles
    for(auto var : triangles) {  
        I = points[var[0]];
        J = points[var[1]];
        L = points[var[2]];
        K = points[var[3]];
      
        right[0] = ((I + J)*0.5)*(J-I);
        mat.setkoeff(1,1,(J-I).koord[0]);
        mat.setkoeff(1,2,(J-I).koord[1]);
        mat.setkoeff(1,3,(J-I).koord[2]);

        right[1] = ((I + L)*0.5)*(L-I);
        mat.setkoeff(2,1,(L-I).koord[0]);
        mat.setkoeff(2,2,(L-I).koord[1]);
        mat.setkoeff(2,3,(L-I).koord[2]);

        right[2] = ((I + K)*0.5)*(K-I);
        mat.setkoeff(3,1,(K-I).koord[0]);
        mat.setkoeff(3,2,(K-I).koord[1]);
        mat.setkoeff(3,3,(K-I).koord[2]);
        
        ergebnis = mat.solve(right);
        p[0] = ergebnis[0];
        p[1] = ergebnis[1];
        p[2] = ergebnis[2];
        circumcenter.push_back(p);
    }

    init_surface();
    ut = kb*T/q;
    ni = sqrt(Nc*Nc*exp((-(Ec-Ev)/ut)));

    std::fill_n(std::back_inserter(volumes), indices.size(), 0);
}


double PSSolver3D::getsurface(long center, long neighbor1) {
    /**
     * \n 
     * -center: first point's index
     * -neighbor: second point's index
     * \n 
     * \n 
     * First the indices of the triangles which contains both points are searched. Then all points circumcenters are ordered so that they build a closed loop.
     * With the class Surface the surface and calculated and then returned.
     */
    std::vector<long> tetraeder, neighbors, copy;
    std::vector<struct List3D> reihenfolge;
    std::vector<long> boundaries;
    std::vector<Vec3d> vertex = {};
    struct List3D new_element; 
    long neighbor_next = -1, neighbor_prev = -1;
    Vec3d Q, N, P1, P2, A, B, AB;
    Surface surf;

   //find all neighbor nodes for center node
    for(unsigned int i = 0; i < triangles.size(); i++) {
        if(std::find(triangles[i].begin(), triangles[i].end(), center) != triangles[i].end()) {
            tetraeder.push_back(i);
            for(unsigned k = 0; k < dim; k++) {
                if(std::find(neighbors.begin(), neighbors.end(), triangles[i][k]) == neighbors.end() && 
                    triangles[i][k] != center) {
                    neighbors.push_back(triangles[i][k]);
                }
           }
       }
    }

    new_element.neighbor_next = -1;
    new_element.neighbor_prev = -1;
    reihenfolge.clear();
    copy = tetraeder;
    
    for (unsigned long z = 0; z < copy.size(); z++) {
        if (std::find(triangles[copy[z]].begin(), triangles[copy[z]].end(), neighbor1) == triangles[copy[z]].end()) {
            copy.erase(copy.begin() + z);
            z = -1;
        }
    }

    //removed all tetraedernumbers in the variable copy, which doesnt include the neighbor[j]
    //index1 = copy.size();
    new_element = circumcenter[copy[0]];
    for(long i = 0; i < 4; i++) {
        if(triangles[copy[0]][i] != center && triangles[copy[0]][i] != neighbor1) {
            if(new_element.neighbor_next == -1) {
                new_element.neighbor_next = triangles[copy[0]][i];
                neighbor_next = triangles[copy[0]][i];
            }
            else {
                new_element.neighbor_prev = triangles[copy[0]][i];
                neighbor_prev = triangles[copy[0]][i];
            }
        }
    }
    reihenfolge.push_back(new_element);
    copy.erase(copy.begin());  

    //hinten dranhänge
    for(unsigned int i = 0; i < copy.size(); i++) {
        for(long t = 0; t < 4; t++) {
            if(triangles[copy[i]][t] == neighbor_next) {
                new_element.neighbor_prev = neighbor_next;
                for(long z = 0; z < 4; z++) {
                    if(triangles[copy[i]][z] != neighbor_next && triangles[copy[i]][z] != center && 
                        triangles[copy[i]][z] != neighbor1) {
                        
                        neighbor_next = triangles[copy[i]][z];
                        new_element.neighbor_next = neighbor_next;
                        break;
                    }
                }
                new_element =  circumcenter[copy[i]];
                reihenfolge.push_back(new_element);
                copy.erase(copy.begin() + i);
                i = -1;
                break;
            }  
        }
    }

    //vorne dranhängen
    std::reverse(std::begin(reihenfolge), std::end(reihenfolge));
    for(unsigned int i = 0; i < copy.size(); i++) {
        for(long t = 0; t < 4; t++) {
            if(triangles[copy[i]][t] == neighbor_prev) {
                new_element.neighbor_next = neighbor_prev;
                for(long z = 0; z < 4; z++) {
                    if(triangles[copy[i]][z] != neighbor_prev && triangles[copy[i]][z] != center &&
                        triangles[copy[i]][z] != neighbor1) {
                        
                        neighbor_prev = triangles[copy[i]][z];
                        new_element.neighbor_prev = neighbor_prev;
                        break;
                    }
                }
                new_element =  circumcenter[copy[i]];
                reihenfolge.push_back(new_element);
                copy.erase(copy.begin() + i);
                i = -1;
                break;
            }  
        }
    }

    std::reverse(std::begin(reihenfolge), std::end(reihenfolge));
    if(copy.size() != 0) std::cout << "fehler size " << copy.size()<< std::endl;   //copy muss size 0 haben

    
    //first and last point are not the same, add buffer points to get a closed loop
    if(reihenfolge[0].neighbor_prev != reihenfolge[reihenfolge.size()-1].neighbor_next) {        
        double lambda;
        //find surfaces that contain center and neighbor point, be careful to not get to much surfaces
        for(unsigned i = 0; i < simple_boundaries.size(); i++) {
            A = points[center];
            B = points[neighbor1];
            AB = (A+B)*0.5;
            
            A = simple_boundaries[i].normalvector;
            B = simple_boundaries[i].surfacepoint;
            if(std::abs(A*(B-AB)) < std::numeric_limits<double>::epsilon()) {
                A = simple_boundaries[i].normalvector;
                unsigned j;
                for(j = 0; j < boundaries.size(); j++) {
                    B = simple_boundaries[boundaries[j]].normalvector;
                    if((A-B).norm() < std::numeric_limits<double>::epsilon()) {
                        break;
                    }
                }
                if(j == boundaries.size()) {
                    boundaries.push_back(i);
                }
            }
        }

        //normal plane
        if(boundaries.size() == 1) {
            if(reihenfolge.size() == 1) std::cout << "Fall vergessen" << std::endl;
            A = simple_boundaries[boundaries[0]].normalvector;
            B = simple_boundaries[boundaries[0]].surfacepoint;
            Q = reihenfolge[reihenfolge.size()-1].p;
            lambda = A*(B-Q);
            B = Q + A*lambda;
            new_element = B;
            reihenfolge.push_back(new_element);

            A = simple_boundaries[boundaries[0]].normalvector;
            B = simple_boundaries[boundaries[0]].surfacepoint;
            Q = reihenfolge[0].p;
            lambda = A*(B-Q);
            B = Q + A*lambda;
            new_element = B;
            std::reverse(std::begin(reihenfolge), std::end(reihenfolge));
            reihenfolge.push_back(new_element);
            std::reverse(std::begin(reihenfolge), std::end(reihenfolge));

            A = points[center];
            B = points[neighbor1];
            A = (A+B)*0.5;
            new_element = A;
            reihenfolge.push_back(new_element); 
            
        }
        else if(boundaries.size() == 2) {
            if(reihenfolge.size() == 1) {
                A = simple_boundaries[boundaries[0]].normalvector;
                B = simple_boundaries[boundaries[0]].surfacepoint;
                Q = reihenfolge[reihenfolge.size()-1].p;
                lambda = A*(B-Q);
                B = Q + A*lambda;
                new_element = B;
                reihenfolge.push_back(new_element);

                A = simple_boundaries[boundaries[1]].normalvector;
                B = simple_boundaries[boundaries[1]].surfacepoint;
                Q = reihenfolge[0].p;
                lambda = A*(B-Q);
                B = Q + A*lambda;
                new_element = B;
                std::reverse(std::begin(reihenfolge), std::end(reihenfolge));
                reihenfolge.push_back(new_element);
                std::reverse(std::begin(reihenfolge), std::end(reihenfolge));

                A = points[center];
                B = points[neighbor1];
                A = (A+B)*0.5;
                new_element = A;
                reihenfolge.push_back(new_element); 
            }
            else if(reihenfolge.size() == 2) {
                A = reihenfolge[0].p;
                B = reihenfolge[0].p;
                if(std::abs((A-B).norm()) < std::numeric_limits<double>::epsilon()) {
                    
                    A = simple_boundaries[boundaries[0]].normalvector;
                    B = simple_boundaries[boundaries[0]].surfacepoint;
                    Q = reihenfolge[reihenfolge.size()-1].p;
                    lambda = A*(B-Q);
                    B = Q + A*lambda;
                    new_element = B;
                    reihenfolge.push_back(new_element);

                    A = simple_boundaries[boundaries[1]].normalvector;
                    B = simple_boundaries[boundaries[1]].surfacepoint;
                    Q = reihenfolge[0].p;
                    lambda = A*(B-Q);
                    B = Q + A*lambda;
                    new_element = B;
                    std::reverse(std::begin(reihenfolge), std::end(reihenfolge));
                    reihenfolge.push_back(new_element);
                    std::reverse(std::begin(reihenfolge), std::end(reihenfolge));

                    A = points[center];
                    B = points[neighbor1];
                    A = (A+B)*0.5;
                    new_element = A;
                    reihenfolge.push_back(new_element); 
                }
                else {
                    std::cout << "Error 3, please contact stella@iue.tuwien.ac.at" << std::endl;
                }
            }
            else {
                std::cout << "Error 1, please contact stella@iue.tuwien.ac.at" << std::endl; 
            }
        }
        else {
            std::cout << "Error 2, please contact stella@iue.tuwien.ac.at" << std::endl; 
        }
    }
        
    // Q = reihenfolge[0].p;
    P2 = points[neighbor1];
    P1 = points[center];
    Q = (P2 + P1)*0.5;  
    N = P2 - P1;
    N = N* (1/N.norm()); 
    for(unsigned int s = 0; s < reihenfolge.size(); s++) {
        P2 = reihenfolge[s].p;
        vertex.push_back(P2);  
    }
    surf.init(reihenfolge.size());
    surf.import_vertex(vertex);
    vertex.clear();
    return surf.get_surface(N);
}
 

double PSSolver3D::getvolume(long center) {
    /**
     * \n
     * center: the index of the point of which the control volume needs to be calculated
     * \n 
     * \n 
     * First all neighbor points are searched. With the function getsurface the surfaces for all neighbor points are calculated.
     */
    std::vector<long> tetraeder, neighbors, copy;
    std::vector<struct List3D> reihenfolge;
    std::vector<Vec3d> vertex = {};
    struct List3D new_element; 
    double  volume = 0;
    Vec3d Q, N, P1, P2;
    Surface surf;
   
    for(unsigned int i = 0; i < triangles.size(); i++) {
       if(std::find(triangles[i].begin(), triangles[i].end(), center) != triangles[i].end()) {
           tetraeder.push_back(i);
           if(std::find(neighbors.begin(), neighbors.end(), triangles[i][0]) == neighbors.end() && 
                triangles[i][0] != center) {
                    neighbors.push_back(triangles[i][0]);
            }
            if(std::find(neighbors.begin(), neighbors.end(), triangles[i][1]) == neighbors.end() && 
                triangles[i][1] != center) {
                    neighbors.push_back(triangles[i][1]);
            }
            if(std::find(neighbors.begin(), neighbors.end(), triangles[i][2]) == neighbors.end() && 
                triangles[i][2] != center) {
                    neighbors.push_back(triangles[i][2]);
            }
            if(std::find(neighbors.begin(), neighbors.end(), triangles[i][3]) == neighbors.end() && 
                triangles[i][3] != center) {
                    neighbors.push_back(triangles[i][3]);
            }
       }
    }
  
    //For all neighbors
    for(unsigned int j = 0; j < neighbors.size(); j++) {
        
       // Q = reihenfolge[0].p;
        P2 = points[neighbors[j]];
        P1 = points[center];
        Q = (P2 + P1)*0.5;  
        N = P2 - P1;
        N = N* (1/N.norm()); 
       // cout << N.norm() << endl;
        for(unsigned int s = 0; s < reihenfolge.size(); s++) {
            P2 = reihenfolge[s].p;
            vertex.push_back(P2);  
        }
        surf.init(reihenfolge.size());
        surf.import_vertex(vertex);
        vertex.clear();
        volume += (Q*N)*getsurface(center, neighbors[j]);
    }
    for(unsigned int s = 0; s < simple_boundaries.size(); s++) {
        if(std::find(simple_boundaries[s].points_indices.begin(), simple_boundaries[s].points_indices.end(), center) != simple_boundaries[s].points_indices.end()) {
            N = simple_boundaries[s].normalvector;
            Q = points[center];
             volume +=(Q*N)*simple_boundaries[s].getvolume(center);
        }
    }

    if(volume < 0) volume = volume * (-1);
    return volume/3;
}


double PSSolver3D::getdistance(long center, long neighbor) {
    /**
     * \n 
     * -center: index of the first point
     * -neighbor: index of the second point
     * \n 
     * \n 
     * The distance between the points with the indices center and neighbor is calculated. The index is referred
     *  to the std::vector points which stores all points
     */
    Vec3d I,J;
    I = points[center];
    J = points[neighbor];
    return std::sqrt((I-J)*(I-J));
}

void PSSolver3D::init_surface() {
    /**
     * \n 
     * Search for all triangles which are on the edges an add them to simple_boundaries. For every boundary plane the normalvector is calculated.
     */
    Vec3d I,J,K, normalvector;
    std::array<double,3> A,B,C;
    std::array<long, 3> reihe;
    long counter1,counter2, counter3, counter4;
    mesh2d mesh; 
    bool status;      
    
    /*
    //check if boundary (only area contains surface) 
        if boundary: calculate normalvector  (3 points are needed even though more exists)
                    add normalvector and a random point on the surface to mesh
                    add this boundary the vector of boundaries (simple_boundaries)
    */
    for(auto it = simple_segments.cbegin(); it != simple_segments.cend(); it++) {
        
        for(unsigned int i = 0; i < it->second.size(); i++) {
            status = true;
            
            for(auto it2 = simple_segments.cbegin(); it2 != simple_segments.cend(); it2++) {
                if(std::find(it2->second.cbegin(), it2->second.cend(), 
                        it->second[i]) != it2->second.cend() && it != it2) {
                    status = false;
                    break;
                }
                
            }
            
            
            if(status == true) {
                I  = polygonpoints[surfacepolygons[it->second[i]][0]];
                J  = polygonpoints[surfacepolygons[it->second[i]][1]];
                K  = polygonpoints[surfacepolygons[it->second[i]][2]];
                normalvector =  (J-I)^(K-J);
                normalvector = normalvector*(1/normalvector.norm());
                mesh.normalvector = normalvector;
                mesh.surfacepoint = J;
                simple_boundaries.push_back(mesh);
            }
            
        }
        
    }
   
    //find point which is not on one of the boundaries
    for(unsigned int i = 0; i < indices.size(); i++) {
        unsigned int j = 0;
        I = points[indices[i]];
        for(; j < simple_boundaries.size(); j++) {
            if(std::abs(I*simple_boundaries[j].normalvector - simple_boundaries[j].normalvector*simple_boundaries[j].surfacepoint) < std::numeric_limits<double>::epsilon()) break;
        }
        if(j == simple_boundaries.size()) break;
    } 
    middlepoint = I;

    //make sure normalvector is outward pointing
    for(unsigned int s = 0; s< simple_boundaries.size(); s++) {
        if(simple_boundaries[s].normalvector * (simple_boundaries[s].surfacepoint - middlepoint) < 0) simple_boundaries[s].normalvector = simple_boundaries[s].normalvector * (-1); 
    }

    /*
        find all triangles which are lying on the boundaries and add them the simple boundaries
        some boundaries planes could be positioned equally. counter1, counter2, counter3, counter4 prevent for 
        adding one triangle to different boundary planes
    */
    for(unsigned int i = 0; i < triangles.size(); i++) {
        counter1 = 0;
        counter2 = 0;
        counter3 = 0;
        counter4 = 0;
        for(unsigned int j = 0; j < simple_boundaries.size(); j++) {
            I = points[triangles[i][0]];
            J = points[triangles[i][1]];
            K = points[triangles[i][2]];
            if( std::abs(I*simple_boundaries[j].normalvector - simple_boundaries[j].normalvector*simple_boundaries[j].surfacepoint) < std::numeric_limits<double>::epsilon() &&
                std::abs(J*simple_boundaries[j].normalvector - simple_boundaries[j].normalvector*simple_boundaries[j].surfacepoint) < std::numeric_limits<double>::epsilon() &&
                std::abs(K*simple_boundaries[j].normalvector - simple_boundaries[j].normalvector*simple_boundaries[j].surfacepoint) < std::numeric_limits<double>::epsilon() && counter1 == 0) {
                
                A[0] = points[triangles[i][0]][0];
                A[1] = points[triangles[i][0]][1];
                A[2] = points[triangles[i][0]][2];
                    
                B[0] = points[triangles[i][1]][0];
                B[1] = points[triangles[i][1]][1];
                B[2] = points[triangles[i][1]][2];
                    
                C[0] = points[triangles[i][2]][0];
                C[1] = points[triangles[i][2]][1];
                C[2] = points[triangles[i][2]][2];

                reihe[0] = triangles[i][0];
                reihe[1] = triangles[i][1];
                reihe[2] = triangles[i][2];
                   
                simple_boundaries[j].init(reihe, A, B, C);
                counter1++;
                    
            }

            I = points[triangles[i][0]];
            J = points[triangles[i][2]];
            K = points[triangles[i][3]];
            if( std::abs(I*simple_boundaries[j].normalvector - simple_boundaries[j].normalvector*simple_boundaries[j].surfacepoint) < std::numeric_limits<double>::epsilon()&&
                std::abs(J*simple_boundaries[j].normalvector - simple_boundaries[j].normalvector*simple_boundaries[j].surfacepoint) < std::numeric_limits<double>::epsilon()&&
                std::abs(K*simple_boundaries[j].normalvector - simple_boundaries[j].normalvector*simple_boundaries[j].surfacepoint) < std::numeric_limits<double>::epsilon() && counter2 == 0) {
    
                A[0] = points[triangles[i][0]][0];
                A[1] = points[triangles[i][0]][1];
                A[2] = points[triangles[i][0]][2];
                    
                B[0] = points[triangles[i][2]][0];
                B[1] = points[triangles[i][2]][1];
                B[2] = points[triangles[i][2]][2];
                    
                C[0] = points[triangles[i][3]][0];
                C[1] = points[triangles[i][3]][1];
                C[2] = points[triangles[i][3]][2];
                   
                reihe[0] = triangles[i][0];
                reihe[1] = triangles[i][2];
                reihe[2] = triangles[i][3];

                simple_boundaries[j].init(reihe, A, B, C);
                counter2++;
            }

            I = points[triangles[i][0]];
            J = points[triangles[i][1]];
            K = points[triangles[i][3]];
               
            if( std::abs(I*simple_boundaries[j].normalvector - simple_boundaries[j].normalvector*simple_boundaries[j].surfacepoint) < std::numeric_limits<double>::epsilon() &&
                std::abs(J*simple_boundaries[j].normalvector - simple_boundaries[j].normalvector*simple_boundaries[j].surfacepoint) < std::numeric_limits<double>::epsilon() &&
                std::abs(K*simple_boundaries[j].normalvector - simple_boundaries[j].normalvector*simple_boundaries[j].surfacepoint) < std::numeric_limits<double>::epsilon() && counter3 == 0) {
                    
                A[0] = points[triangles[i][0]][0];
                A[1] = points[triangles[i][0]][1];
                A[2] = points[triangles[i][0]][2];
                    
                B[0] = points[triangles[i][1]][0];
                B[1] = points[triangles[i][1]][1];
                B[2] = points[triangles[i][1]][2];
                    
                C[0] = points[triangles[i][3]][0];
                C[1] = points[triangles[i][3]][1];
                C[2] = points[triangles[i][3]][2];

                reihe[0] = triangles[i][0];
                reihe[1] = triangles[i][1];
                reihe[2] = triangles[i][3];
                    
                simple_boundaries[j].init(reihe, A, B, C);
                counter3++;
            }

            I = points[triangles[i][1]];
            J = points[triangles[i][2]];
            K = points[triangles[i][3]];
                
            if( std::abs(I*simple_boundaries[j].normalvector - simple_boundaries[j].normalvector*simple_boundaries[j].surfacepoint) < std::numeric_limits<double>::epsilon() &&
                std::abs(J*simple_boundaries[j].normalvector - simple_boundaries[j].normalvector*simple_boundaries[j].surfacepoint) < std::numeric_limits<double>::epsilon() &&
                std::abs(K*simple_boundaries[j].normalvector - simple_boundaries[j].normalvector*simple_boundaries[j].surfacepoint) < std::numeric_limits<double>::epsilon() && counter4 == 0) {
                    
                A[0] = points[triangles[i][1]][0];
                A[1] = points[triangles[i][1]][1];
                A[2] = points[triangles[i][1]][2];
                    
                B[0] = points[triangles[i][2]][0];
                B[1] = points[triangles[i][2]][1];
                B[2] = points[triangles[i][2]][2];
                    
                C[0] = points[triangles[i][3]][0];
                C[1] = points[triangles[i][3]][1];
                C[2] = points[triangles[i][3]][2];
                    
                reihe[0] = triangles[i][1];
                reihe[1] = triangles[i][2];
                reihe[2] = triangles[i][3];

                simple_boundaries[j].init(reihe, A, B, C);
                counter4++;
            }
        }        
    }
}

void PSSolver3D::set_polygonpoints(std::vector<std::array<double,3> > points) { 
    polygonpoints = points;
}

void PSSolver3D::set_surfacepolygon(std::vector<std::vector<long> >  points) {
    surfacepolygons = points;
}

void PSSolver3D::test() {
    double volume = 0;
    create_mesh();
    for(unsigned i = 0; i < indices.size(); i++) volume += getvolume(indices[i]);
    std::cout << "Gesamtvolumen " << volume << std::endl;
}

void PSSolver3D::add_simple_boundary(std::vector<long> boundary, std::string name) {
    simple_segments.insert(std::make_pair(name, boundary));
}


void PSSolver3D::scalepoints(std::vector<std::vector<double>> &grid, double factor) {
    for(auto it = grid.begin(); it !=grid.end(); it++) {
        for(auto it2 = it->begin(); it2 != it->end(); it2++) {
            *it2 = *it2*factor;
        }
    }

    for(auto it = polygonpoints.begin(); it != polygonpoints.end(); it++) {
        for(auto it2 = it->begin(); it2 != it->end(); it2++) {
            *it2 = *it2*factor;
        }
    }
}
