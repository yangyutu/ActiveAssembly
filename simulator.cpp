#include"simulator.h"


Simulator::Simulator(std::shared_ptr<Model> model0,std::shared_ptr<Controller> controller0):
                    model(model0),controller(controller0){
    controlFrequency = 1.0/model->dt();
    assignmentFrequency = 10;
    nstep = 1000;
}

void Simulator::run(){
    
    model->createInitialState();
    controller->calAssignment(model->getCurrState());
    controller->calControl(model->getCurrState());
    int iter;
    for(int s=0; s < nstep; s++){
        if ((s+1)%assignmentFrequency == 0){
            controller->calAssignment(model->getCurrState());
        }
        controller->calControl(model->getCurrState());        
        model->run(controlFrequency);
    }
}