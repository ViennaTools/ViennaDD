#include "functions.hpp"


void Surface::init(int length) { 
  numbers = length; 
}

void Surface::import_vertex(std::vector<Vec3d> data) { 
  vertex = data; 
}

double Surface::get_surface(Vec3d n) {
  /**
   * \n 
   * -n: normalvector of the plane
   * \n 
   * \n 
   * The points need to be ordered to calculate the right surface. 
   */
  double surface = 0;
  Vec3d  result;
  result.koord[0] = 0;
  result.koord[1] = 0;
  result.koord[2] = 0;
  for (long j = 0; j < numbers; j++) {
    if (j == numbers-1) {
       result = result + (vertex[j] ^ vertex[0]);
    }
     
    else {
      result = result + (vertex[j] ^ vertex[j + 1]);
    }
      
  }
  surface = result * n;
  if (surface < 0)
    surface *= -1;
  return surface / 2;
}


Vec2d Vec2d::operator*(double scalar) {
  Vec2d result;
  result.koord[0] = koord[0] * scalar;
  result.koord[1] = koord[1] * scalar;
  return result;
}

Vec2d Vec2d::operator+(Vec2d p2) {
  Vec2d result;
  result.koord[0] = koord[0] + p2.koord[0];
  result.koord[1] = koord[1] + p2.koord[1];
  return result;
}

Vec2d Vec2d::operator-(Vec2d p2) {
  Vec2d result;
  result.koord[0] = koord[0] - p2.koord[0];
  result.koord[1] = koord[1] - p2.koord[1];
  return result;
}

void Vec2d::operator=(std::array<double, 2> p2) {
  koord = p2;
}

void Vec2d::operator=(std::vector<double> p2) {
  koord[0] = p2[0];
  koord[1] = p2[1];
}
    
double Vec2d::operator*(Vec2d p2) {
  return koord[0] * p2.koord[0] +  koord[1] * p2.koord[1];
}

Vec3d Vec3d::operator*(double scalar) {
  Vec3d result;
  result.koord[0] = koord[0] * scalar;
  result.koord[1] = koord[1] * scalar;
  result.koord[2] = koord[2] * scalar;
  return result;
}

Vec3d Vec3d::operator+(Vec3d p2) {
  Vec3d result;
  result.koord[0] = koord[0] + p2.koord[0];
  result.koord[1] = koord[1] + p2.koord[1];
  result.koord[2] = koord[2] + p2.koord[2];
  return result;
}

Vec3d Vec3d::operator-(Vec3d p2) {
  Vec3d result;
  result.koord[0] = koord[0] - p2.koord[0];
  result.koord[1] = koord[1] - p2.koord[1];
  result.koord[2] = koord[2] - p2.koord[2];
  return result;
}

void Vec3d::operator=(std::array<double, 3> p2) {
  koord = p2;
}

void Vec3d::operator=(std::vector<double> p2) {
  koord[0] = p2[0];
  koord[1] = p2[1];
  koord[2] = p2[2];
}
    
double Vec3d::operator*(Vec3d p2) {
  return koord[0] * p2.koord[0] +  koord[1] * p2.koord[1] +  koord[2] * p2.koord[2];
}

double Vec3d::norm() {
  return std::sqrt(koord[0]*koord[0] + koord[1]*koord[1] + koord[2]*koord[2]);
}


Vec3d Vec3d::operator^(Vec3d B) {
	Vec3d res;
	res.koord[0] = koord[1]*B.koord[2] - koord[2]*B.koord[1];
	res.koord[1] = koord[2]*B.koord[0] - koord[0]*B.koord[2];
	res.koord[2] = koord[0]* (B.koord[1]) - koord[1]* (B.koord[0]);
	return res;
}




void Matrix2d::setkoeff(long row, long column, double value) {
  koeff[(row-1)*2 + column -1] = value;
}

std::array<double,2> Matrix2d::solve(std::array<double, 2> vec) {
  double det = 0;
  double buffer;
  std::array<double, 2> result;
  det = koeff[0]*koeff[3] - koeff[1] * koeff[2];
  
  //det C1
  buffer = vec[0]*koeff[3] - koeff[1] * vec[1];
  result[0] = buffer/det;
        
  //det C2
  buffer = koeff[0]*vec[1] - vec[0] * koeff[2];
  result[1] = buffer/det;
  return result;
}


void Matrix3d::setkoeff(long row, long column, double value) {
  koeff[(row-1)*3 + column -1] = value;
}


std::array<double,3> Matrix3d::solve(std::array<double, 3> vec) {
  double det = 0;
  double buffer;
  std::array<double, 3> result;
         
  det =   koeff[0]*koeff[4]*koeff[8] +
          koeff[1]*koeff[5]*koeff[6] + 
          koeff[2]*koeff[3]*koeff[7] - 
          koeff[2]* koeff[4] * koeff[6] - 
          koeff[1]* koeff[3] * koeff[8] - 
          koeff[0]* koeff[5] * koeff[7];
                                       
  //det C1
  buffer =  vec[0]*koeff[4]*koeff[8] +
            koeff[1]*koeff[5]*vec[2] + 
            koeff[2]*vec[1]*koeff[7] - 
            koeff[2]* koeff[4] * vec[2] - 
            koeff[1]* vec[1] * koeff[8] - 
            vec[0]* koeff[5] * koeff[7];
  result[0] = buffer/det;

  //det C2
  buffer =  koeff[0]*vec[1]*koeff[8] +
            vec[0]*koeff[5]*koeff[6] + 
            koeff[2]*koeff[3]*vec[2] - 
            koeff[2]* vec[1] * koeff[6] - 
            vec[0]* koeff[3] * koeff[8] - 
            koeff[0]* koeff[5] * vec[2];
  result[1] = buffer/det;

  //det C3
  buffer =  koeff[0]*koeff[4]*vec[2] +
            koeff[1]*vec[1]*koeff[6] + 
            vec[0]*koeff[3]*koeff[7] - 
            vec[0]* koeff[4] * koeff[6] - 
            koeff[1]* koeff[3] * vec[2] - 
            koeff[0]* vec[1] * koeff[7];
  result[2] = buffer/det;
        
  return result;
}


void mesh2d::init(std::array<long, 3> tr, std::array<double,3> p1,
          std::array<double,3> p2, std::array<double,3> p3) {
  /**
   * \n 
   * -tr: contain the indices of the 3 points whcih are connected to one triangle
   * -p1: coordinates of the first point
   * -p2: coordinates of the second point
   * -p3: coordinates of the third point
   * \n 
   * \n 
   * At first new points are stored and then the circumcenter of the new triangle is calculated
   */
  Vec3d I, L, J;
  Matrix3d mat;
  std::array <double,3> right, ergebnis;
  for(long j = 0; j < 3; j++) {
    if(std::find(points_indices.begin(), points_indices.end(), tr[j]) == points_indices.end()) {
                //not saved
      points_indices.push_back(tr[j]);
      if(j == 0) points.push_back(p1);
      else if(j == 1)points.push_back(p2);
      else points.push_back(p3);
    }
  }
  
  triangles.push_back(tr);
        
  I = p1;
  J = p2;
  L = p3;
 

  //Calculate circumcenter
  right[2] = normalvector*J;
  mat.setkoeff(3,1,normalvector.koord[0]);
  mat.setkoeff(3,2,normalvector.koord[1]);
  mat.setkoeff(3,3,normalvector.koord[2]);

  right[0] = ((I + J)*0.5)*(J-I);
  mat.setkoeff(1,1,(J-I).koord[0]);
  mat.setkoeff(1,2,(J-I).koord[1]);
  mat.setkoeff(1,3,(J-I).koord[2]);

  right[1] = ((L + I)*0.5)*(L-I);
  mat.setkoeff(2,1,(L-I).koord[0]);
  mat.setkoeff(2,2,(L-I).koord[1]);
  mat.setkoeff(2,3,(L-I).koord[2]);

  ergebnis = mat.solve(right);
  circumcenter.push_back(ergebnis);
}

double mesh2d::getvolume(long center) {
  /**
   * \n 
   * -center: the index of which the controlvolume needs to be calculated
   * \n 
   * \n 
   * First alle triangles which contain the center are searched. Then the circumcenters are stored in the right. The class Surface 
   * is used to calculate the surface the control volume
   */
  std::vector<long> neighbors = {}, fig = {}, copy = {};
  std::vector<struct List3D> reihenfolge;
  struct List3D new_element;
  std::array<double, 3> new_point;
  std::vector<Vec3d> vertex = {};
  Vec3d A,B, AB, L;
  Surface surf;
  long neighbor_next = -1, neighbor_prev = -1, buffer;
  
  //find all triangles that contains center node
  for(unsigned long j = 0; j < triangles.size(); j++) {
    if(triangles[j][0] == center || triangles[j][1] == center || triangles[j][2] == center) {
      fig.push_back(j);
    }
  }
  
  //add first point and find both neighbor nodes
  new_element.p = circumcenter[fig[0]];
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
  //delete first triangle and add element  
  reihenfolge.push_back(new_element);
  fig.erase(fig.begin());
        

  //add new elements at the end
  for(unsigned long i = 0; i < fig.size(); i++) {
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
        new_element.p = circumcenter[fig[i]];
        reihenfolge.push_back(new_element);
        fig.erase(fig.begin() + i);
        i = -1;
        break;
      }  
    }
  }


  //add new elements at the front    
  std::reverse(std::begin(reihenfolge), std::end(reihenfolge));
  for(unsigned long i = 0; i < fig.size(); i++) {
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
        new_element.p = circumcenter[fig[i]];
        reihenfolge.push_back(new_element);
        fig.erase(fig.begin() + i);
        i = -1;
        break;
      }  
    }
  }
  std::reverse(std::begin(reihenfolge), std::end(reihenfolge));
   
  //if not closed loop add buffer node to close the loop (vertice)
  if(reihenfolge[0].neighbor_prev != reihenfolge[reihenfolge.size()-1].neighbor_next) {
    buffer = std::distance(points_indices.begin(), std::find(points_indices.begin(), points_indices.end(), center));
    A = points[buffer];
    buffer = std::distance(points_indices.begin(), std::find(points_indices.begin(), points_indices.end(), reihenfolge[reihenfolge.size()-1].neighbor_next));
    B = points[buffer];
    AB = (A + B) * 0.5;
    new_element = AB;
    reihenfolge.push_back(new_element);


    buffer = std::distance(points_indices.begin(), std::find(points_indices.begin(), points_indices.end(), center));
    A = points[buffer];
    buffer = std::distance(points_indices.begin(), std::find(points_indices.begin(), points_indices.end(), reihenfolge[0].neighbor_prev));
    B = points[buffer];
    AB = (A + B) * 0.5;
    new_element = AB;

    std::reverse(std::begin(reihenfolge), std::end(reihenfolge));
    reihenfolge.push_back(new_element);
    std::reverse(std::begin(reihenfolge), std::end(reihenfolge));

    buffer = std::distance(points_indices.begin(), std::find(points_indices.begin(), points_indices.end(), center));
    new_point[1] = points[buffer][1];
    new_point[0] = points[buffer][0];
    new_point[2] = points[buffer][2];
    new_element.p = new_point;
    reihenfolge.push_back(new_element);  
  }
    
  //add all points to vertex and then calculate surface
  for(unsigned long s = 0; s < reihenfolge.size(); s++) {
    L = reihenfolge[s].p;
    vertex.push_back(L);  
  }

  surf.init(reihenfolge.size());
  surf.import_vertex(vertex);
  vertex.clear();

  return surf.get_surface(normalvector);
}
