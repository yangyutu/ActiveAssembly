#pragma once
#include<vector>
#include<memory>
#include<random>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
//#include "CellList.h"

class CellList; 
typedef std::shared_ptr<CellList> cellList_ptr;
class Model {
public:

    struct particle {
        double r[3],F[3],ori_vec[3][3];
        double phi;
        double theta;
        double u;
        int targetIdx;
        double cost;
    };
    typedef std::shared_ptr<particle> particle_ptr;
    typedef std::vector<particle_ptr> state;
   
    Model(){}
    Model(int np, int dim0,double radius0, std::string filetag0, std::string ini, std::string targetFile, cellList_ptr cell);
    ~Model() {
    }
    void run();
    void run(int steps);
    void createInitialState();
    state getCurrState(){return particles;}
    int getDimP(){return dimP;}
    std::vector<Model::particle> getTargets(){return targets;}
    double dt(){return dt_;}
    int np(){return numP;}
private:
    void calForces();
    void calForcesHelper(int i, int j, double F[3]);
    bool cellListFlag;
    std::shared_ptr<CellList> cellList;
    int dimP;
    static const double kb, T, vis;
    int numP;
    double radius;
    double LJ,rm;
    double eps[3][3][3];
    std::vector<double> velocity={0.0,2.0e-6,5.0e-6};
    std::vector<particle_ptr> particles;
    std::vector<int> control;
    std::string iniFile;
    double dt_, cutoff, mobility, diffusivity_r, diffusivity_t;
    std::default_random_engine rand_generator;
    std::shared_ptr<std::normal_distribution<double>> rand_normal;
    int trajOutputInterval;
    int timeCounter,fileCounter;
    std::ofstream trajOs, opOs;
    std::string filetag;
    std::vector<Model::particle> targets;
    void outputTrajectory(std::ostream& os);
    void outputOrderParameter(std::ostream& os);
    void readxyz(const std::string filename);
    void updateBodyFrameVec();
    void readTarget(std::string filename);
    void getPermutator();
};



