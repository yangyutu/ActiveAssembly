#include<iostream>
#include "simulator.h"
#include "controller.h"
#include "model.h"
#include "CellList.h"

int main(){

    int N = 79;
    int dim = 2;
    double radius = 1e-6;
    cellList_ptr cell(new CellList(3.0*radius,3,20,300.0*radius,300.0*radius,300.0*radius));
    std::shared_ptr<Model> m(new Model(N,dim,radius,"control","config2d.txt" ,"triangle_79p.txt",nullptr));
    std::shared_ptr<Controller> c(new Controller(N,dim, radius, m->getTargets()));
    Simulator simulator(m,c);
    simulator.run();  
    return 0;
}

