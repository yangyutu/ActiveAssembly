#include<iostream>
#include "simulator.h"
#include "controller.h"
#include "model.h"
#include "CellList.h"

int main(){

    int N = 2022;
    int dim = 3;
    double radius = 1e-6;
    cellList_ptr cell(new CellList(3.0*radius,2,10,300.0*radius,300.0*radius,300.0*radius));
    std::shared_ptr<Model> m(new Model(N,dim,radius,"control","config3d.txt" ,"toporing_big.txt",cell));
    std::shared_ptr<Controller> c(new Controller(N,dim, radius, m->getTargets()));
    Simulator simulator(m,c);
    simulator.run();  
    return 0;
}

