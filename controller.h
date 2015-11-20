#pragma once
#include<vector>
#include<map>
#include "model.h"

class Controller {
public:

    Controller(int dim, double radius);

    ~Controller() {
    }
    typedef std::vector<int> control;
    void calControl(Model::state s);
    void calAssignment(Model::state s);
    void getErrorDist();

private:
    int dim;
    double radius;
    void readErrorMap();
    void readPolicyMap();
    void readTargets();
    std::vector<Model::particle2d> targets;
    std::vector<int> assignment;
};
