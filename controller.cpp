#include<cmath>
#include "controller.h"
#include "HungarianAlg.h"
#include<fstream>
Controller::Controller(int dim0, double radius0):dim(dim0){
//    readErrorMap();
//    initialize targets
    radius = radius0;
    double R = (int)(dim/2.0/M_PI) * radius;
    double d_theta = 2.0*M_PI/dim;
    R= 10;
    for(int i = 0; i < dim; i++){
        targets.push_back(Model::particle2d());
        targets[i].r[0] = R*cos(d_theta*i);
        targets[i].r[1] = R*sin(d_theta*i);
    }
    


    
//    for(int i = 0; i < dim; i++){
//        targets[i].r[0] *= radius;
 //       targets[i].r[1] *= radius;
 //   }    
    std::ifstream is;
    is.open("triangle_79p.txt");
    std::string line;
    std::stringstream linestream;
    double dum;
    for (int i = 0; i < dim; i++) {
        getline(is, line);
        std::stringstream linestream(line);
        linestream >> dum;
        linestream >> targets[i].r[0];
        linestream >> targets[i].r[1];
    }
    
    
    
    std::ofstream os;
    os.open("target.txt");
    for(int i = 0; i < dim; i++){
        os << i << "\t";
        os << targets[i].r[0] << "\t";
        os << targets[i].r[1] << std::endl;
    }    
}

void Controller::calControl(Model::state s){
//    Controller::control control;
    for(int i=0; i < dim; i++){
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

void Controller::calAssignment(Model::state s) {
    vector< vector<double> > Cost(dim, vector<double>(dim));
    for(int i=0; i<dim; i++){
	for(int j=0; j<dim; j++){
            double rx = targets[j].r[0] - s[i]->r[0]/radius ;
            double ry = targets[j].r[1] - s[i]->r[1]/radius ;
            double proj_x = cos(s[i]->phi)*rx + sin(s[i]->phi)*ry;
            double proj_y = -sin(s[i]->phi)*rx + cos(s[i]->phi)*ry;
            double c = pow((proj_x - 2)/1.2,2.0) + pow((proj_y),2.0);
            Cost[i][j] = c;
        }
    }   
    AssignmentProblemSolver APS;
    APS.Solve(Cost, assignment);
    for(int i=0; i < dim; i++){
        s[i]->targetIdx = assignment[i];
        s[i]->cost = Cost[i][assignment[i]];
    }
}