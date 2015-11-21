#include"simulator.h"


Simulator::Simulator(std::shared_ptr<Model> model0,std::shared_ptr<Controller> controller0):
                    model(model0),controller(controller0){
    controlFrequency = 1.0/model->dt();
    assignmentFrequency = 10;
    nstep = 1000;
}

void Simulator::run(){
    double totalCost;   
    model->createInitialState();
    totalCost = controller->calAssignment(model->getCurrState(),model->getDimP());
    controller->calControl(model->getCurrState(),model->getDimP());
    int iter;
    std::cout << totalCost << std::endl;
    for(int s=0; s < nstep; s++){
        if ((s+1)%assignmentFrequency == 0){
            totalCost = controller->calAssignment(model->getCurrState(),model->getDimP());
            std::cout << totalCost << std::endl;
        }
        controller->calControl(model->getCurrState(),model->getDimP());        
        model->run(controlFrequency);
    }
}

