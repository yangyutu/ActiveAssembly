#include<iostream>
#include "simulator.h"
#include "controller.h"
#include "model.h"


int main(){

    int N = 79;
    double radius = 1e-6;
    std::shared_ptr<Model> m(new Model(N,radius,"control"));
    std::shared_ptr<Controller> c(new Controller(N, radius));
    Simulator simulator(m,c);
    simulator.run();  
    return 0;
}
