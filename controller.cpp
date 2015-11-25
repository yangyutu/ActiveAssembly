#include<cmath>
#include "controller.h"
#include "HungarianAlg.h"
#include<fstream>
Controller::Controller(int numP0,int dim0, double radius0,
        std::vector<Model::particle> target0):numP(numP0), dimP(dim0), targets(target0){
//    readErrorMap();
//    initialize targets
    radius = radius0;
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

void Controller::calControl(Model::state s, int dim){
    if (dim ==2){
        calControl2d(s);    
    } else if(dim ==3){
        calControl3d(s);
    }
}


double Controller::calAssignment(Model::state s, int dim){
    double totalCost;
    if (dim ==2){
      totalCost = calAssignment2d(s);    
    } else if(dim ==3){
       totalCost = calAssignment3d(s);
    }
    return totalCost;
}


void Controller::calControl2d(Model::state s){
//    Controller::control control;
    for(int i=0; i < numP; i++){
        int t_idx = assignment[i];
        double rx = targets[t_idx].r[0] - s[i]->r[0]/radius ;
        double ry = targets[t_idx].r[1] - s[i]->r[1]/radius ;
        double dot_prod = cos(s[i]->phi)*rx + sin(s[i]->phi)*ry;
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


double Controller::calAssignment3d(Model::state s) {
    vector< vector<long> > Cost(numP, vector<long>(numP));
    double totalCost = 0.0;
    for(int i=0; i<numP; i++){
	for(int j=0; j<numP; j++){
            double rx = targets[j].r[0] - s[i]->r[0]/radius ;
            double ry = targets[j].r[1] - s[i]->r[1]/radius ;
            double rz = targets[j].r[2] - s[i]->r[2]/radius;
            double proj_x = s[i]->ori_vec[0][0]*rx + s[i]->ori_vec[0][0]*ry+s[i]->ori_vec[0][0]*rz;
            double proj_y = s[i]->ori_vec[0][1]*rx + s[i]->ori_vec[0][1]*ry+s[i]->ori_vec[0][1]*rz;
            double proj_z = s[i]->ori_vec[0][2]*rx + s[i]->ori_vec[0][2]*ry+s[i]->ori_vec[0][2]*rz;;
            double c = pow((proj_x - 2)/1.2,2.0) + pow((proj_y),2.0) + pow((proj_z),2.0);
            Cost[i][j] = (long)(sqrt(c)*10.0);
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

void Controller::calControl3d(Model::state s){
//    Controller::control control;
    for(int i=0; i < numP; i++){
        int t_idx = assignment[i];
        double rx = targets[t_idx].r[0] - s[i]->r[0]/radius ;
        double ry = targets[t_idx].r[1] - s[i]->r[1]/radius ;
        double rz = targets[t_idx].r[2] - s[i]->r[2]/radius;
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


double Controller::calAssignment2d(Model::state s) {
    vector< vector<long> > Cost(numP, vector<long>(numP));
    double totalCost = 0.0;
    for(int i=0; i<numP; i++){
	for(int j=0; j<numP; j++){
            double rx = targets[j].r[0] - s[i]->r[0]/radius ;
            double ry = targets[j].r[1] - s[i]->r[1]/radius ;
            double proj_x = cos(s[i]->phi)*rx + sin(s[i]->phi)*ry;
            double proj_y = -sin(s[i]->phi)*rx + cos(s[i]->phi)*ry;
            double c = pow((proj_x - 2)/1.2,2.0) + pow((proj_y),2.0);
            Cost[i][j] = (long)(sqrt(c)*10.0);
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