#include<cmath>
#include<set>
#include<limits>
#include "controller.h"
#include "HungarianAlg.h"
#include<fstream>
#include "common.h"
#include <lemon/dijkstra.h>

//#include "icpPointToPoint.h"
//#include "icpPointToPlane.h"
#include <armadillo>

using namespace lemon;
extern Parameter parameter;

Controller::Controller(Model::state targets){
//    readErrorMap();
//    initialize targets
    targets_ = targets;
    radius = parameter.radius;
    numP = parameter.N;
    dimP = parameter.dim;
    
    colliInfo.vSet[0] = 0.0;
    colliInfo.vSet[1] = 2.0;
    colliInfo.vSet[2] = 5.0;
    colliInfo.colliThresh = 6.0;
    for (int i = 0; i < numP; i++)
        availControl.push_back(2);
    
   
    min_x = -25;
    min_y = -25;
    del_x = 0.5;
    del_y = 0.5;
    x_binNum = 101;
    y_binNum = 101;
    
    
    
     std::array<Array2D_type::index, 2> dim = {x_binNum,y_binNum};
    maps[0] = std::make_shared<Array2D_type>(dim);
    maps[1] = std::make_shared<Array2D_type>(dim);
    maps[2] = std::make_shared<Array2D_type>(dim);
    
    this->readCostMap();
    
    if (parameter.landmarkFlag){
        this->constructLandmark();
    
    }
            /*
    double R = (int)(dim/2.0/M_PI) * radius;
    double d_theta = 2.0*M_PI/dim;
    R= 10;
    for(int i = 0; i < dim; i++){
        targets.push_back(Model::particle());
        targets[i].r[0] = R*cos(d_theta*i);
        targets[i].r[1] = R*sin(d_theta*i);
    }
    
*/
}

void Controller::calControl(Model::state s, Model::state targets, int dim){
    if (dim ==2){
        calControl2d(s, targets);    
    } else if(dim ==3){
        calControl3d(s, targets);
    }
}


double Controller::calAssignment(Model::state s, Model::state targets, int dim){
    double totalCost;
    if (dim ==2){
      totalCost = calAssignment2d(s, targets);    
      if (parameter.motionFlag){
//          totalCost = calAssignment2d
      }
    } else if(dim ==3){
       totalCost = calAssignment3d(s, targets);
    }
    return totalCost;
}


void Controller::calControl2d(Model::state s, Model::state targets){
//    Controller::control control;
    for(int i=0; i < numP; i++){
        int t_idx = s[i]->targetIdx;
        if (t_idx >=0){
            double rx = targets[t_idx]->r[0] - s[i]->r[0]/radius ;
            double ry = targets[t_idx]->r[1] - s[i]->r[1]/radius ;
            double dot_prod = cos(s[i]->phi)*rx + sin(s[i]->phi)*ry;
            if (dot_prod < 0) {
                s[i]->u = 0;
            } else {
                if (dot_prod < 1) {
                    s[i]->u = 0;
                    ;
                } else if (dot_prod < 3.7) {
                    s[i]->u = 1;
                } else {
                    s[i]->u = 2;
                }
            }

            if (s[i]->u > availControl[i]) {
                s[i]->u = availControl[i];
            }
        }
    }
}

void Controller::readCostMap(){    
    std::ifstream is;
    std::string line;
    for (int mapIdx = 0; mapIdx < 3; mapIdx++) {
        std::stringstream ss;
        ss << mapIdx;
        is.open("costMap"+ss.str()+".txt");
        for (int i = 0; i < x_binNum; i++) {
            getline(is, line);
            std::stringstream linestream(line);
            for (int j = 0; j < y_binNum; j++){
                linestream >> (*(maps[mapIdx]))[i][j];
//                std::cout<<(*(maps[mapIdx]))[i][j]<<"\t";
            }
                
        }
        is.close();
    }
}

double Controller::getCostFromMap(double x,double y, double z, int mapIndex){
    int z_bin, x_bin, y_bin;
//
//    if (dim == 2) {
//        z_bin = 0;
//        x_bin = (int) ((x - min_x) / del_x);
//        y_bin = (int) ((y - min_y) / del_y);
//    } else {
        x_bin = (int) ((x - min_x) / del_x);
        y_bin = (int) ((y - min_y) / del_y);
//        z_bin = (int) ((z - min_z) / del_z);
       
        double c = (*(maps[mapIndex]))[x_bin][y_bin];
        return c;
}

/*
void Controller::calAvoidance2d(Model::state s){

    for (int i=0; i < numP; i++){
        availControl[i] = 2;
    }
    
    for (int i=0; i < numP-1; i++){
        for (int j = i+1; j < numP; j++){
            double v_dot = cos(s[i]->phi)*cos(s[j]->phi)+
            sin(s[i]->phi)*sin(s[j]->phi);
            if (v_dot > sqrt(3.0)/2.0) break;  // if in the same direction within 30 degree difference, no need to avoid          
//          otherwise calculate the intersection
            double r[3];
            double v[3];
            r[0] = (s[i]->r[0] - s[j]->r[0])/radius;
            r[1] = (s[i]->r[1] - s[j]->r[1])/radius;
            r[2] = (s[i]->r[2] - s[j]->r[2])/radius;
            double dist = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
            dist = sqrt(dist);             
            if (dist > colliInfo.colliThresh) break; // too far away to hit
            
            v[0] = cos(s[i]->phi) - cos(s[j]->phi);
            v[1] = sin(s[i]->phi) - sin(s[j]->phi);
            v[2] = 0;
            double b = r[0] * v[0] + r[1] * v[1] + r[2] * v[2];
            if (b > 0) break;// collision cannot happen since particles are moving away from each other
                       
            double InvertMat[2][2];
            double det = -cos(s[i]->phi)*sin(s[j]->phi)+cos(s[j]->phi)*sin(s[i]->phi);
            InvertMat[0][0] = -sin(s[j]->phi);
            InvertMat[1][1] = cos(s[i]->phi);
            InvertMat[0][1] = cos(s[j]->phi);
            InvertMat[1][0] = -sin(s[i]->phi);
            
            double t0,t1;
            if (det < 1e-6){ // this is the case two particle are almost anti-parallel
                double perpDist = abs(-r[0]*sin(s[i]->phi) + r[1]*cos(s[j]->phi));
                if (perpDist > 2.0) break; // the perpendicular distance is big, they won't collide
                t0 = 0.5*dist;
                t1 = 0.5*dist;
            } else {
                t0 = (InvertMat[0][0]*b[0]+InvertMat[0][1]*b[1])/det;
                t1 = (InvertMat[1][0]*b[0]+InvertMat[1][1]*b[1])/det;
            }
            // if both greater than zero then they intersect on their front path
            if (t1 > 0 && t0 > 0){
                if (t1 > t0){
                    
                    while(colliInfo.vSet[availControl[j]] >= t1){
                        availControl[j]--;
                        if (availControl[j] == 0) break;
                    }
                } else {
                    while(colliInfo.vSet[availControl[i]] >= t0){
                        availControl[i]--;
                        if (availControl[i] == 0) break;
                    }
                }           
            } else {// t1*t0<0, where one hit another in the back, then reduce the speed of particle that hits others
                if (t1 > 0){
                    
                    while(colliInfo.vSet[availControl[j]] >= t1){
                        availControl[j]--;
                        if (availControl[j] == 0) break;
                    }
                } else if(t0 >0) {
                    while(colliInfo.vSet[availControl[i]] >= t0){
                        availControl[i]--;
                        if (availControl[i] == 0) break;
                    }
                }            
            }
 
        }
    }
                
}
*/
void Controller::calAvoidance2d_simpleCollision(Model::state s){
    for (int i=0; i < numP; i++){
        availControl[i] = 2;
    }
    
    for (int i=0; i < numP; i++){
        for (int j = 0; j < numP; j++){
            //if one is diffuse only
//            std::cout << j << "\t";
            if (i != j){
            if (availControl[i] == 0) break;
            double r[3];
            double v[3];
            r[0] = (s[i]->r[0] - s[j]->r[0])/radius;
            r[1] = (s[i]->r[1] - s[j]->r[1])/radius;
            r[2] = (s[i]->r[2] - s[j]->r[2])/radius;
            double dist = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
            dist = sqrt(dist);             
            if (dist > colliInfo.colliThresh) continue;
            int totalControl = availControl[i];
            int iter = 0;
            while (iter < totalControl){
                iter++;
                // assume the second particle is static
                v[0] = colliInfo.vSet[availControl[i]]*cos(s[i]->phi) - 0.0;
                v[1] = colliInfo.vSet[availControl[i]]*sin(s[i]->phi) - 0.0;
                v[2] = 0;
                double v_sq = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
                double b = r[0] * v[0] + r[1] * v[1] + r[2] * v[2];
                if (b > 0) break;// collision cannot happen at this speed
                double delta = b * b - v_sq * (dist * dist - 4.0);
                if (delta < 0) break; // collision cannot happen at this speed
                // if collision can happen
                double t = (-b - sqrt(delta)) / v_sq;
                if (t < 1.0) {
                    availControl[i]--;   
                } else {
                    break;
                }
            }
            }
        }
    }
    for (int i=0; i < numP; i++){
        s[i]->availControl=availControl[i];
//        std::cout << i << "\t" << s[i]->availControl << std::endl;
    }
}

/*

void Controller::calAvoidance2d_colliVersion(Model::state s){
// this version use collision checking in Molecular simulation
// I have not figured it out yet
    for (int i=0; i < numP; i++){
        availControl[i] = 2;
    }
    
    for (int i=0; i < numP-1; i++){
        for (int j = i+1; j < numP; j++){
            double v_dot = cos(s[i]->phi)*cos(s[j]->phi)+
            sin(s[i]->phi)*sin(s[j]->phi);
            if (v_dot > 0) break;  // if in the same direction, no need to avoid          
            if (availControl[i] == 0 && availControl[j] == 0) break;
            double r[3];
            double v[3];
            r[0] = (s[i]->r[0] - s[j]->r[0])/radius;
            r[1] = (s[i]->r[1] - s[j]->r[1])/radius;
            r[2] = (s[i]->r[2] - s[j]->r[2])/radius;
            double dist = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
            dist = sqrt(dist);             
            if (dist > colliInfo.colliThresh) break;
            int totalControl = availControl[i] + availControl[j];
            int iter = 0;
            while (iter < totalControl){
                iter++;
                v[0] = colliInfo.vSet[availControl[i]]*cos(s[i]->phi) - colliInfo.vSet[availControl[i]]*cos(s[j]->phi);
                v[1] = colliInfo.vSet[availControl[i]]*sin(s[i]->phi) - colliInfo.vSet[availControl[j]]*sin(s[j]->phi);
                v[2] = 0;
                double v_sq = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
                double b = r[0] * v[0] + r[1] * v[1] + r[2] * v[2];
                if (b > 0) break;// collision cannot happen at this speed
                double delta = b * b - v_sq * (dist * dist - 4.0);
                if (delta < 0) break; // collision cannot happen at this speed
                // if collision can happen
                double t = (-b - sqrt(delta)) / v_sq;
                if (t < 1.0) {
                    if (availControl[i] = availControl[j]) {
                        int randomSelect;
                        randomSelect = rand() % 2;
                        if (randomSelect == 0) {
                            availControl[i]--;
                        } else {
                            availControl[j]--;
                        }
                    } else if (availControl[i] > availControl[j]) {
                        availControl[i]--;
                    } else {
                        availControl[j]--;
                    }
                }
            }
        }
    }
                
}
  */  
double Controller::calAssignment3d(Model::state s, Model::state targets) {
    vector< vector<double> > Cost(numP, vector<double>(numP));
    double totalCost = 0.0;
    for(int i=0; i<numP; i++){
	for(int j=0; j<numP; j++){
            double rx = targets[j]->r[0] - s[i]->r[0]/radius ;
            double ry = targets[j]->r[1] - s[i]->r[1]/radius ;
            double rz = targets[j]->r[2] - s[i]->r[2]/radius;
            double proj_x = s[i]->ori_vec[0][0]*rx + s[i]->ori_vec[0][0]*ry+s[i]->ori_vec[0][0]*rz;
            double proj_y = s[i]->ori_vec[0][1]*rx + s[i]->ori_vec[0][1]*ry+s[i]->ori_vec[0][1]*rz;
            double proj_z = s[i]->ori_vec[0][2]*rx + s[i]->ori_vec[0][2]*ry+s[i]->ori_vec[0][2]*rz;
            double c = pow((proj_x - 2)/1.2,2.0) + pow((proj_y),2.0) + pow((proj_z),2.0);
//            Cost[i][j] = (long)(sqrt(c)*10.0);
            Cost[i][j] = c;
        }
    }   
    AssignmentProblemSolver APS;
    APS.Solve(Cost, assignment);
    for(int i=0; i < numP; i++){
        s[i]->targetIdx = assignment[i];
        s[i]->cost = Cost[i][assignment[i]];
        totalCost += s[i]->cost;
    }
    return totalCost;
}

void Controller::calControl3d(Model::state s, Model::state targets){
//    Controller::control control;
    for(int i=0; i < numP; i++){
        int t_idx = assignment[i];
        double rx = targets[t_idx]->r[0] - s[i]->r[0]/radius ;
        double ry = targets[t_idx]->r[1] - s[i]->r[1]/radius ;
        double rz = targets[t_idx]->r[2] - s[i]->r[2]/radius;
        double dot_prod = s[i]->ori_vec[0][0]*rx + s[i]->ori_vec[0][1]*ry+s[i]->ori_vec[0][2]*rz;
        if (dot_prod < 0){
            s[i]->u = 0;
        } else {
            if (dot_prod < 1){
                s[i]->u = 0;;
            } else if(dot_prod < 2.8){
                s[i]->u = 1;
            } else{
                s[i]->u = 2;
            }
        }
    }
}


double Controller::calAssignment2d(Model::state s, Model::state targets) {
    this->calAvoidance2d_simpleCollision(s);
    
    vector< vector<double> > Cost(numP, vector<double>(numP));
    double totalCost = 0.0;
    for(int i=0; i<numP; i++){
	for(int j=0; j<numP; j++){
            double rx = targets[j]->r[0] - s[i]->r[0]/radius ;
            double ry = targets[j]->r[1] - s[i]->r[1]/radius ;
            double proj_x = cos(s[i]->phi)*rx + sin(s[i]->phi)*ry;
            double proj_y = -sin(s[i]->phi)*rx + cos(s[i]->phi)*ry;
            
            double dist = proj_x*proj_x + proj_y*proj_y;
            dist = sqrt(dist);
            double c;
            if (dist > 20){
                c = pow((proj_x - 2)/1.2,2.0) + pow((proj_y),2.0);
            } else{
                c = this->getCostFromMap(proj_x,proj_y,0,availControl[i]);
            }
                Cost[i][j] = c;
//            Cost[i][j] = (long)(sqrt(c)*10.0);
        }
    }   
    AssignmentProblemSolver APS;
    APS.Solve(Cost, assignment);
    for(int i=0; i < numP; i++){
        s[i]->targetIdx = assignment[i];
        s[i]->cost = Cost[i][assignment[i]];
        targets[s[i]->targetIdx]->targetIdx = i;
        totalCost += s[i]->cost;
    }
    return totalCost;
}

void Controller::translate_2d(double phi,Model::state s){
// first select a subset of particles that is in the right direction
//  let them go until shape change too much
    this->calInlier(s);
    for (int i = 0; i < numP; i++){
        double proj;
        proj = cos(s[i]->phi - phi);
        if (proj > sqrt(3.0)/2.0 && s[i]->nbcount >=6) {
            s[i]->u = 2;
        }
    }
    
}

void Controller::rotate_2d(Model::state s){
// first determine a rotational center
    double center_s[3];
    this->calInlier(s);    
    this->calWeightCenter(s,center_s,1);
    for (int k = 0; k < 3; k++){
        std::cout << "rotate center_s: " << center_s[k]/radius << std::endl;
    }
    

    for (int i = 0; i < numP; i++){
        double proj;
        double pos_phi = atan2(s[i]->r[1]/radius-center_s[1]/radius,s[i]->r[0]/radius-center_s[0]/radius) + M_PI/2.0;
        proj = cos(s[i]->phi - pos_phi);
        if (proj > sqrt(3.0)/2.0 && s[i]->nbcount >=6) {
            s[i]->u = 2;
        }
    }

    
}

void Controller::calWeightCenter(Model::state s, double center[3],int flag) {
    if (flag == 0) {
        for (int k = 0; k < 3; k++) {
            center[k] = 0;
            for (int i = 0; i < numP; i++) {
                center[k] += s[i]->r[k];
            }
            center[k] /= numP;
        }
    } else {
        for (int k = 0; k < 3; k++) {
            center[k] = 0;
            double weight_sum = 0.0;
            for (int i = 0; i < numP; i++) {
                double weight = s[i]->inlier;
                center[k] += s[i]->r[k] * weight;
                weight_sum += weight;
            }
            center[k] /= (weight_sum);
        }
    }
}

void Controller::alignTarget_t(Model::state s, Model::state targets){

    double center_target[3], center_s[3];
        
    this->calWeightCenter(targets,center_target,0);
    for (int k = 0; k < 3; k++){
        
        std::cout << "center_target: " << center_target[k] << std::endl;
    }

    this->calWeightCenter(s,center_s,1);
    for (int k = 0; k < 3; k++){
        center_s[k] /= (radius);
        std::cout << "center_s: " << center_s[k] << std::endl;
        
    }

    // for do the shift
    for (int i = 0; i < numP; i++) {
        for (int k = 0; k < 3; k++) {
            targets[i]->r[k] += center_s[k] - center_target[k];
        }
    }
}

/*
void Controller::register_2d(Model::state s, Model::state targets){
  // define a 3 dim problem with 10000 model points
  // and 10000 template points:
  int32_t dim = 2;
  int32_t num = numP;

  // allocate model and template memory
  double* M = (double*)calloc(dim*num,sizeof(double));
  double* T = (double*)calloc(dim*num,sizeof(double));

  // set model and template points
  for (int i = 0; i < numP; i++){
      T[i*dim+0] = targets[i]->r[0];
      T[i*dim+1] = targets[i]->r[1];
  
      M[i*dim+0] = s[i]->r[0]/radius;
      M[i*dim+1] = s[i]->r[1]/radius;
  }
  // start with identity as initial transformation
  // in practice you might want to use some kind of prediction here
  Matrix R = Matrix::eye(2);
  Matrix t(2,1);
  
  
      double center_target[3], center_s[3];
        
    this->calWeightCenter(targets,center_target,0);
    for (int k = 0; k < 3; k++){
        
        std::cout << "center_target: " << center_target[k] << std::endl;
    }

    this->calWeightCenter(s,center_s,1);
    for (int k = 0; k < 3; k++){
        center_s[k] /= (radius);
        std::cout << "center_s: " << center_s[k] << std::endl;
        
    }

    // for do the shift
    for (int i = 0; i < numP; i++) {
        for (int k = 0; k < 3; k++) {
            t(0,0)= -center_s[0] + center_target[0];
            t(1,0)= -center_s[1] + center_target[1];
        }
    }
  

  // run point-to-plane ICP (-1 = no outlier threshold)
  cout << endl << "Running ICP (point-to-point, no outliers)" << endl;
  IcpPointToPoint icp(M,num,dim);
  icp.fit(T,num,R,t,2.5);

  // results
  cout << endl << "Transformation results:" << endl;
  cout << "R:" << endl << R << endl << endl;
  cout << "t:" << endl << t << endl << endl;

  // free memory
  free(M);
  free(T);

  for (int i = 0; i < numP; i++){
      targets[i]->r[0] -= t(0,0) ;
      targets[i]->r[1] -= t(1,0);
      targets[i]->r[0] = R(0,0)*targets[i]->r[0] + R(0,1)*targets[i]->r[1];
      targets[i]->r[1] = R(1,0)*targets[i]->r[0] + R(1,1)*targets[i]->r[1];      
  }
  
  
  
}
*/
void Controller::calInlier(Model::state s){

    double r[3],dist;
    int size = s.size();
    for (int i = 0; i < size; i++) {
        s[i]->nbcount = 0;
        s[i]->inlier = 0;
    }
    for (int i = 0; i < size - 1; i++) {
        for (int j = i + 1; j < size; j++) {
            dist = 0.0;
            for (int k = 0; k < dimP; k++) {

                r[k] = (s[j]->r[k] - s[i]->r[k]) / radius;
                dist += pow(r[k], 2.0);
            }
            dist = sqrt(dist);
            if (dist < 2.5) {
                s[i]->nbcount++;
                s[j]->nbcount++;
            }

        }
    }
    for (int i = 0; i < size; i++) {
       if( s[i]->nbcount >=3){
           s[i]->inlier = 1;
       }
    }
    
}

double Controller::calSeqAssignment(Model::state s, Model::state targets, int expandFlag){

    if (expandFlag) {
        deviation = 0.0;
        double r[3];
        for (int i = 0; i < numP; i++) {
            if (targets[i]->marked) {

                r[0] = targets[i]->r[0] - s[targets[i]->targetIdx]->r[0] / radius;
                r[1] = targets[i]->r[1] - s[targets[i]->targetIdx]->r[1] / radius;
                r[2] = targets[i]->r[2] - s[targets[i]->targetIdx]->r[2] / radius;
                deviation += sqrt(pow(r[0], 2) + pow(r[1], 2) + pow(r[2], 2));
            }
            
        }
        deviation /= numMarked;

        if (deviation < 1.0) {
            this->expandTargets();
        }
    }
    this->calAvoidance2d_simpleCollision(s);
    
    vector< vector<double> > Cost(numP, vector<double>(numMarked));

    for(int i=0; i<numP; i++){
	for(int j=0; j<numMarked; j++){
            double rx = targets[markedIdx[j]]->r[0] - s[i]->r[0]/radius ;
            double ry = targets[markedIdx[j]]->r[1] - s[i]->r[1]/radius ;
            double proj_x = cos(s[i]->phi)*rx + sin(s[i]->phi)*ry;
            double proj_y = -sin(s[i]->phi)*rx + cos(s[i]->phi)*ry;
            
            double dist = proj_x*proj_x + proj_y*proj_y;
            dist = sqrt(dist);
            double c;
            if (dist > 20){
                c = pow((proj_x - 2)/1.2,2.0) + pow((proj_y),2.0);
            } else{
                c = this->getCostFromMap(proj_x,proj_y,0,availControl[i]);
            }
                Cost[i][j] = c;
//            Cost[i][j] = (long)(sqrt(c)*10.0);
        }
    }   
    AssignmentProblemSolver APS;
    APS.Solve(Cost, assignment);
    double totalCost = 0.0;
    for(int i=0; i < numP; i++){
        if(assignment[i] >= 0){
            s[i]->targetIdx = markedIdx[assignment[i]];
            s[i]->cost = Cost[i][assignment[i]];
            targets[markedIdx[assignment[i]]]->targetIdx = i;
            totalCost += s[i]->cost;
        } else {
            availControl[i] =0; // no need to control;
        
        }
    }
    totalCost /=numMarked;
    return totalCost;
}

// here is alignement with correponding
// specific derivation can be found at thesis: Evaluation of surface registration algorithms for PET motion correction

void Controller::alignTarget_rt(Model::state s, Model::state targets){

    int size = s.size();
    
    this->calInlier(s);
//    NRmatrix<double> b(2, 2);
    arma::mat cov(2,2);
    double center_s[3], center_target[3];
    
    this->calWeightCenter(targets,center_target,0);
    for (int k = 0; k < 3; k++){
        
        std::cout << "center_target: " << center_target[k] << std::endl;
    }

    this->calWeightCenter(s,center_s,1);
    for (int k = 0; k < 3; k++){
        center_s[k] /= (radius);
        std::cout << "center_s: " << center_s[k] << std::endl;
        
    }
    
    
    for (int i =0; i <2; i++){
        for (int j = 0; j < 2; j++) {
            cov(i,j) = 0;
            for (int k = 0; k < size; k++){
                int t_idx = s[k]->targetIdx;
               double weight = s[i]->inlier;
                cov(i,j) += (s[k]->r[j]/radius - center_s[j])*(targets[t_idx]->r[i] - center_target[i])*weight;            
            }
        }
    }
    arma::mat U;
    arma::vec sin;
    arma::mat V;

    arma::svd(U,sin,V,cov);
    
    arma::mat R;
    R = V*U.t();
    R.print("rotation matrix:");

    double tran[2];
    tran[0] = center_s[0] - (R(0,0)*center_target[0]+R(0,1)*center_target[1]);
    tran[1] = center_s[1] - (R(1,0)*center_target[0]+R(1,1)*center_target[1]);   
    std::cout << "translation:" << std::endl;
    std::cout << tran[0] << "\t" << tran[1] << std::endl;
    for (int i = 0; i < size; i++) {
        targets[i]->r[0] = R(0,0)*targets[i]->r[0] + R(0,1)*targets[i]->r[1];
        targets[i]->r[1] = R(1,0)*targets[i]->r[0] + R(1,1)*targets[i]->r[1];

        targets[i]->r[0] += tran[0];
        targets[i]->r[1] += tran[1];        
    }
    
}

void Controller::translateCargo_2d(double phi, Model::state s){
    this->calInlier(s);
    for (int i = 0; i < numP; i++){
        double proj;
        proj = cos(s[i]->phi - phi);
        if (proj > sqrt(3.0)/2.0 && s[i]->nbcount >=4) {
            s[i]->u = 2;
        }
    }
    // the cargo just diffuse
    s[0]->u = 0;
}
void Controller::alignCargo(Model::state s,Model::state targets){
    double t[3];
  //  s[0] is the cargo, targets[0] is the center of the template
    t[0] = s[0]->r[0]/radius - targets[0]->r[0];
    t[1] = s[0]->r[1]/radius - targets[0]->r[1];
    t[2] = s[0]->r[2]/radius - targets[0]->r[2];
    
    std::cout << "align with cargo " << std::endl;
    std::cout << t[0] << "\t" << t[1] << std::endl;
    
    for (int i = 0; i < numP; i++){
        targets[i]->r[0] += t[0];
        targets[i]->r[1] += t[1]; 
        targets[i]->r[2] += t[2]; 
    
    }
    
}

void Controller::buildTargetGraph(){
    numMarked = 0;
    
    marked_t = std::make_shared<lemon::ListGraph::NodeMap<int>>(targetG);
    index_t = std::make_shared<lemon::ListGraph::NodeMap<int>>(targetG);
    for (int i = 0; i < numP; i++){
        nodes_t.push_back(targetG.addNode());
        (*marked_t)[nodes_t[i]] = targets_[i]->marked;
       (*index_t)[nodes_t[i]] = i;
       if (targets_[i]->marked){
        numMarked++;
       markedIdx.push_back(i);        
       }

    }
    numSurface=numMarked;
// add edges
   double r[3],dist;

    for (int i = 0; i < numP - 1; i++) {
        for (int j = i + 1; j < numP; j++) {
            dist = 0.0;
            for (int k = 0; k < dimP; k++) {

                r[k] = (targets_[j]->r[k] - targets_[i]->r[k]);
                dist += pow(r[k], 2.0);
            }
            dist = sqrt(dist);
            if (dist < 2.5) {
                edges_t.push_back(targetG.addEdge(nodes_t[i],nodes_t[j]));
            }

        }
    }
   std::cout << "target graph info:" << std::endl;
   std::cout << "nodes: " << countNodes(targetG) << std::endl;
   std::cout << "nodes: " << countArcs(targetG) << std::endl;
   
}

void Controller::expandTargets(){
    
    std::set<ListGraph::Node> nodes;
    for (int i = 0; i < numP; i++){
        if ((*marked_t)[nodes_t[i]]){
            for (ListGraph::OutArcIt e(targetG,nodes_t[i]); e != INVALID;++e){
                ListGraph::Node n = targetG.oppositeNode(nodes_t[i],e);
                if ((*marked_t)[n] == 0){
                    nodes.insert(n);
                }
 
            }
        }
    }
    
    for (ListGraph::Node n: nodes){
        numMarked++;
        (*marked_t)[n] = 1;
        int idx = (*index_t)[n];
        targets_[idx]->marked = 1; 
        markedIdx.push_back(idx); 
    }
//    numSurface = nodes.size();
}

void Controller::constructLandmark() {
    landmarkLength = parameter.landmarkLength;
    landmarkDist = parameter.landmarkDist;
    
    landmark_pos = std::make_shared<lemon::ListGraph::NodeMap<Model::particle_ptr>>(landmarkG);
    internalLength = std::make_shared<lemon::ListGraph::EdgeMap<double>>(landmarkG);
    length = std::make_shared<lemon::ListGraph::EdgeMap<double>>(landmarkG);
    int count = 0;
    for (int i = 0; i < landmarkLength; i++) {
        for (int j = 0; j < landmarkLength; j++) {
            nodes_l.push_back(landmarkG.addNode());
            Model::particle_ptr p =  std::make_shared<Model::particle>(i*landmarkDist - 0.5 * landmarkLength,
                    j*landmarkDist - 0.5 * landmarkLength, 0.0);
            (*landmark_pos)[nodes_l[count]] = p;
            count++;
        }
    }
    
    numLandmark = count + 1;
    for (int i = 0; i < numLandmark; i++) {
        shortPathMat.push_back(std::vector<double>(numLandmark,0.0));
    }
    for (int i = 0; i < numLandmark - 1; i++) {
        for (int j = i + 1; j < numLandmark; j++) {
            double dx = (*landmark_pos)[nodes_l[i]]->r[0] - (*landmark_pos)[nodes_l[j]]->r[0];
            double dy = (*landmark_pos)[nodes_l[i]]->r[1] - (*landmark_pos)[nodes_l[j]]->r[1];
            double d = sqrt(dx * dx + dy * dy);
            if (d < sqrt(3) * landmarkDist) {
                edges_l.push_back(landmarkG.addEdge(nodes_l[i], nodes_l[j]));
                (*internalLength)[edges_l[edges_l.size() - 1]] = d;
                (*length)[edges_l[edges_l.size() - 1]] = 0.0;
            }
        }
    }
    
    
    
    
}

void Controller::assignLandmarkIdx(Model::state s, double scale) {
    for (int i = 0; i < numP; i++) {
        s[i]->nbLandmark.clear();
        s[i]->nbLandmarkDist.clear();
        for (int j = 0; j < numLandmark; j++) {
            double dx = (*landmark_pos)[nodes_l[j]]->r[0] - s[i]->r[0] / scale;
            double dy = (*landmark_pos)[nodes_l[j]]->r[1] - s[i]->r[1] / scale;
            double d = sqrt(dx * dx + dy * dy);
            if (d < sqrt(2) * landmarkDist) {
                s[i]->nbLandmark.push_back(j);
                s[i]->nbLandmarkDist.push_back(d);
            }
        }
    }
}

double Controller::calExtraCost(Model::state s, double r1[3], double scale1, 
                            double r2[3], double scale2){
    
    r1[0] *= scale1;
    r1[1] *= scale1;
    r1[2] *= scale1;
    
    r2[0] *= scale2;
    r2[1] *= scale2;
    r2[2] *= scale2;
    
    
    double vec1[3];
    double vec2[3];
    
    vec1[0] = r2[0] - r1[0];
    vec1[1] = r2[1] - r1[1];
    vec1[2] = r2[2] - r1[2];
    double norm1 = sqrt(vec2[0]*vec2[0] + vec2[1]*vec2[1] + vec2[2]*vec2[2]);
    vec1[0] /= norm1;
    vec1[1] /= norm1;
    vec1[2] /= norm1;
    
    for (int i = 0; i < numP; i++){
        vec2[0] = s[i]->r[0] - r1[0];
        vec2[1] = s[i]->r[1] - r1[1];
        vec2[2] = s[i]->r[2] - r1[0];
        
        double proj = vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
        double norm = sqrt(vec2[0]*vec2[0] + vec2[1]*vec2[1] + vec2[2]*vec2[2]);
        double dist = sqrt(norm*norm - proj);
        
        if (dist < (2.0*radius) && proj > radius && proj < (norm-radius)){
            return blockCost;
        }
        
           
    }

    return 0.0;
}

void Controller::calShortestPathMat(Model::state s, Model::state targets) {
    
    this->assignLandmarkIdx(s, radius);
    this->assignLandmarkIdx(targets,1.0);
    
    for (int i = 0; i < numLandmark; i++){
        for (ListGraph::OutArcIt e(landmarkG,nodes_l[i]); e != INVALID;++e){
            ListGraph::Node n = landmarkG.oppositeNode(nodes_l[i],e);
            double extraCost = calExtraCost(s,(*landmark_pos)[nodes_l[i]]->r,radius,(*landmark_pos)[n]->r,radius);
            (*length)[e] = (*internalLength)[e] + extraCost;
        }
    }
    
    for (int i = 0; i < numP; i++){
        for (int j=0; s[i]->nbLandmark.size(); j++){
            int idx = s[i]->nbLandmark[j];
            double extraCost = calExtraCost(s,s[i]->r,1.0,(*landmark_pos)[nodes_l[idx]]->r,radius);
            s[i]->nbLandmarkDist[j] += extraCost;
        }
    }
    
    for (int i = 0; i < numP; i++){
        for (int j: targets[i]->nbLandmark){
            int idx = s[i]->nbLandmark[j];
            double extraCost = calExtraCost(s,targets[i]->r,radius,(*landmark_pos)[nodes_l[j]]->r,radius);
            targets[i]->nbLandmarkDist[j] += extraCost;
        }
    }
    
 
    Dijkstra<ListGraph, ListGraph::EdgeMap<double>> dij(landmarkG, *length);
    for (int i = 0; i < numLandmark; i++) {
        dij.run(nodes_l[i]);
        for (int j = 0; j < numLandmark; j++) {
            shortPathMat[i][j] = dij.dist(nodes_l[j]);
        }
    }
}

double Controller::calPathDistViaLandmark(Model::state s, int i, Model::state targets, int j){

    double r[3];
    r[0] = s[i]->r[0]/radius - targets[i]->r[0];
    r[1] = s[i]->r[1]/radius - targets[i]->r[1];
    r[2] = s[i]->r[2]/radius - targets[i]->r[2];
    
    // calculate the shortpathmat
    this->calShortestPathMat(s,targets);
        
    double pathDist = sqrt(pow(r[0],2)+pow(r[1],2)+pow(r[2],2)) + 
        this->calExtraCost(s,s[i]->r,1.0,targets[i]->r,radius);

    for (int ii = 0; ii < s[i]->nbLandmark.size(); ii++){
        for (int jj=0; jj < targets[j]->nbLandmark.size();jj++){
            int idx1 = s[i]->nbLandmark[ii];
            int idx2 = s[j]->nbLandmark[jj];
            double pathDistTemp = s[i]->nbLandmarkDist[ii] + targets[j]->nbLandmarkDist[jj]+
                    shortPathMat[idx1][idx2];
            if (pathDistTemp < pathDist){
                pathDist = pathDistTemp;
                s[i]->landmarkIdx = idx1;
                targets[j]->landmarkIdx = idx2;
            }
            
        }
    }
}
