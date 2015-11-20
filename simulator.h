#include<memory>
#include "model.h"
#include "controller.h"

class Simulator{
public:
    Simulator(){}
    Simulator(std::shared_ptr<Model> model0,std::shared_ptr<Controller> controller0);
    ~Simulator(){}
    void run();
    
private:
    std::shared_ptr<Model> model;
    std::shared_ptr<Controller> controller;
    int nstep;
    int assignmentFrequency;
    int controlFrequency;
};