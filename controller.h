#pragma once
#include<vector>
#include<map>
#include "model.h"

class Controller {
public:

    Controller(int nump, int dimP, double radius, std::vector<Model::particle> targets);

    ~Controller() {
    }
    typedef std::vector<int> control;
    double calAssignment(Model::state s, int dimP);
    void calControl(Model::state s, int dimP);
    void calControl2d(Model::state s);
    double calAssignment2d(Model::state s);
    void calControl3d(Model::state s);
    double calAssignment3d(Model::state s);
    void getErrorDist();

private:
    int dimP, numP;
    double radius;
    void readErrorMap();
    void readPolicyMap();
    void readTargets();
    std::vector<Model::particle> targets;
    std::vector<int> assignment;
};
