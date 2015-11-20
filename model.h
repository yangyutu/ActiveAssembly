#pragma once
#include<vector>
#include<memory>
#include<random>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
class Model {
public:

    struct particle2d {
        double r[2],F[2];
        double phi;
        double u;
        int targetIdx;
        double cost;
    };
    typedef std::shared_ptr<particle2d> particle2d_ptr;
    typedef std::vector<particle2d_ptr> state;

    Model(){}
    Model(int np, double radius0, std::string filetag0);
    ~Model() {
    }
    void run();
    void run(int steps);
    void createInitialState();
    state getCurrState(){return particles;}
    double dt(){return dt_;}
    int np(){return numP;}
private:
    void calForces();
    static const int dimP=2;
    static const double kb, T, vis;
    int numP;
    double radius;
    double LJ,rm;
    std::vector<double> velocity={0.0,2.0e-6,5.0e-6};
    std::vector<particle2d_ptr> particles;
    std::vector<int> control;
    double dt_, cutoff, mobility, diffusivity_r, diffusivity_t;
    std::default_random_engine rand_generator;
    std::shared_ptr<std::normal_distribution<double>> rand_normal;
    int trajOutputInterval;
    int timeCounter,fileCounter;
    std::ofstream trajOs, opOs;
    std::string filetag;
    void outputTrajectory(std::ostream& os);
    void outputOrderParameter(std::ostream& os);
    void readxyz(const std::string filename);
};



