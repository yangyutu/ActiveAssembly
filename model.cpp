#include "model.h"

double const Model::T = 293.0;
double const Model::kb = 1.38e-23;
double const Model::vis = 1e-3;

Model::Model(int np, double radius0, std::string filetag0):numP(np), 
        radius(radius0), filetag(filetag0){
    rand_normal = std::make_shared<std::normal_distribution<double>>(0.0, 1.0);

    for(int i = 0; i < numP; i++){
        particles.push_back(particle2d_ptr(new Model::particle2d));    
    }
    dt_ = 0.001;
    diffusivity_t = 2.145e-13;
    diffusivity_r = 0.2145;
    mobility = diffusivity_t/kb/T;
    trajOutputInterval = 1000;
    fileCounter = 0;
    cutoff = 5.0;
    LJ = 3.0*Model::kb*Model::T/radius;
    rm = 2.3;
            
    
}

void Model::run() {
    if (this->timeCounter == 0 || ((this->timeCounter + 1) % trajOutputInterval == 0)) {
        this->outputTrajectory(this->trajOs);
//        this->outputOrderParameter(this->opOs);
    }
    
    calForces();
//    particles[0]->phi = 0.0;
//    particles[1]->phi = -M_PI;
//    particles[0]->u = 2;
//    particles[1]->u = 2;
    
    
    for (int i = 0; i < numP; i++) {
        
        particles[i]->r[0] += mobility * particles[i]->F[0] * dt_ +
                    velocity[particles[i]->u] * cos(particles[i]->phi) * dt_;
                +   sqrt(2.0 * diffusivity_t * dt_) * (*rand_normal)(rand_generator);
        particles[i]->r[1] += mobility * particles[i]->F[1] * dt_ +
                    velocity[particles[i]->u] * sin(particles[i]->phi) * dt_; 
                +   sqrt(2.0 * diffusivity_t * dt_) * (*rand_normal)(rand_generator);
            
            
        particles[i]->phi += sqrt(2.0 * diffusivity_r * dt_) * (*rand_normal)(rand_generator);
    }
    this->timeCounter++;
    
}

void Model::run(int steps){
    for (int i = 0; i < steps; i++){
	run();
    }
}
void Model::calForces(){
    double r[dimP], dist;
    for (int i = 0; i < numP; i++){
        for (int k = 0; k < dimP; k++){
        particles[i]->F[k] = 0.0;
        }
    }
    
    for(int i = 0; i < numP-1; i++){
        for(int j = i+1; j < numP; j++){
            dist = 0.0;
            for(int k = 0; k < dimP; k++){
                r[k] = (particles[j]->r[k] - particles[i]->r[k])/radius;
                dist += pow(r[k],2.0);
            }
            dist = sqrt(dist);
            if(dist < 2.0){
                std::cout << "overlap " << i << "\t" << j << std::endl;
                dist = 2.0;
            }           
            if (dist < cutoff) {
                double Fpp = LJ*(-12.0*pow((rm/dist),12)/dist+12.0*pow((rm/dist),7)/dist);
                Fpp += -9e-13*exp(-33.3*(dist-2.0));
            for(int k = 0; k < dimP; k++){
                particles[i]->F[k] += Fpp * r[k]/dist;
                particles[j]->F[k] += -Fpp * r[k]/dist;
            }   
            }

        
        }
    }
}

void Model::createInitialState(){

    this->readxyz("config.txt");
    std::stringstream ss;
    std::cout << "model initialize at round " << fileCounter << std::endl;
    ss << this->fileCounter++;
    if (trajOs.is_open()) trajOs.close();
//    if (opOs.is_open()) opOs.close();

    this->trajOs.open(filetag + "xyz_" + ss.str() + ".txt");
//    this->opOs.open(filetag + "op" + ss.str() + ".dat");
    this->timeCounter = 0;

}

void Model::outputTrajectory(std::ostream& os) {

    for (int i = 0; i < numP; i++) {
        os << i << "\t";
        os << particles[i]->r[0]/radius << "\t";
        os << particles[i]->r[1]/radius << "\t";
        os << particles[i]->phi<< "\t";
        os << particles[i]->cost<< "\t";
        os << particles[i]->u<< "\t";
        os << particles[i]->targetIdx<< "\t";
        os << std::endl;
    }
}

void Model::readxyz(const std::string filename) {
    std::ifstream is;
    is.open(filename.c_str());
    std::string line;
    double dum;
    for (int i = 0; i < numP; i++) {
        getline(is, line);
        std::stringstream linestream(line);
        linestream >> dum;
        linestream >> particles[i]->r[0];
        linestream >> particles[i]->r[1];
        linestream >> particles[i]->phi;
    }
    for (int i = 0; i < numP; i++) {
        particles[i]->r[0] *=radius;
        particles[i]->r[1] *=radius;
    }
    
    is.close();
}