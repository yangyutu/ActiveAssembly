#include<iostream>
#include "simulator.h"
#include "controller.h"
#include "model.h"


int main(){

    int N = 782;
    int dim = 3;
    double radius = 1e-6;
    std::shared_ptr<Model> m(new Model(N,dim,radius,"control", "target_topo.txt"));
    std::shared_ptr<Controller> c(new Controller(N,dim, radius, m->getTargets()));
    Simulator simulator(m,c);
    simulator.run();  
    return 0;
}
