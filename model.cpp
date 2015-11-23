#include "model.h"
#include "CellList.h"

double const Model::T = 293.0;
double const Model::kb = 1.38e-23;
double const Model::vis = 1e-3;

Model::Model(int np, int dim0, double radius0, std::string filetag0,std::string ini, std::string targetFile, cellList_ptr cell):numP(np), 
        dimP(dim0), radius(radius0), filetag(filetag0), cellList(cell), iniFile(ini){
    rand_normal = std::make_shared<std::normal_distribution<double>>(0.0, 1.0);

    for(int i = 0; i < numP; i++){
        particles.push_back(particle_ptr(new Model::particle));    
    }
    dt_ = 0.001;
    diffusivity_t = 2.145e-13;
//    diffusivity_t = qe-18;
    diffusivity_r = 0.2145;
    mobility = diffusivity_t/kb/T;
    trajOutputInterval = 1000;
    fileCounter = 0;
    cutoff = 5.0;
    LJ = 3.0*Model::kb*Model::T/radius;
    rm = 2.3;
    for(int i = 0; i < numP; i++){
        targets.push_back(Model::particle());
    }        
    this->readTarget(targetFile);
    
    this->getPermutator();
    
    cellListFlag = false;
    if (cell!=nullptr){
        cellListFlag = true;
    }
    
}

void Model::run() {
    if (this->timeCounter == 0 || ((this->timeCounter + 1) % trajOutputInterval == 0)) {
        this->outputTrajectory(this->trajOs);
//        this->outputOrderParameter(this->opOs);
    }
    if (cellListFlag){
            double builtCount = cellList->buildList(particles);
            if (builtCount!=numP){
                std::cout << "build imcomplete" << std::endl;
            }
    }
    calForces();
//    particles[0]->phi = 0.0;
//    particles[1]->phi = -M_PI;
//    particles[0]->u = 2;
//    particles[1]->u = 2;
    
    if (dimP == 2){
        for (int i = 0; i < numP; i++) {

            particles[i]->r[0] += mobility * particles[i]->F[0] * dt_ +
                        velocity[particles[i]->u] * cos(particles[i]->phi) * dt_;
                    +   sqrt(2.0 * diffusivity_t * dt_) * (*rand_normal)(rand_generator);
            particles[i]->r[1] += mobility * particles[i]->F[1] * dt_ +
                        velocity[particles[i]->u] * sin(particles[i]->phi) * dt_; 
                    +   sqrt(2.0 * diffusivity_t * dt_) * (*rand_normal)(rand_generator);


            particles[i]->phi += sqrt(2.0 * diffusivity_r * dt_) * (*rand_normal)(rand_generator);
        }
    } else if(dimP == 3){
        for (int i = 0; i < numP; i++) {
            for(int j = 0; j < dimP; j++){
                particles[i]->r[j] += mobility * particles[i]->F[j] * dt_ +
                        velocity[particles[i]->u] * particles[i]->ori_vec[0][j]* dt_;
                    +   sqrt(2.0 * diffusivity_t * dt_) * (*rand_normal)(rand_generator);
            }
        // update the orientation vector via spherical surface diffusion
            double phitemp = sqrt(2.0 * diffusivity_r * dt_) * (*rand_normal)(rand_generator);
            double thetatemp = sqrt(2.0 * diffusivity_r * dt_) * (*rand_normal)(rand_generator);
            for(int j = 0; j < dimP; j++){
                particles[i]->ori_vec[0][j] += particles[i]->ori_vec[1][j]*phitemp + 
                    particles[i]->ori_vec[2][j]*thetatemp;
            }
        }
        this->updateBodyFrameVec();  
        
    }
    this->timeCounter++;
    
}

void Model::run(int steps){
    for (int i = 0; i < steps; i++){
	run();
    }
}

void Model::readTarget(std::string filename){
    std::ifstream is;
    is.open(filename);
    std::string line;
    std::stringstream linestream;
    double dum;
    for (int i = 0; i < numP; i++) {
        getline(is, line);
        std::stringstream linestream(line);
        linestream >> dum;
        linestream >> targets[i].r[0];
        linestream >> targets[i].r[1];
        linestream >> targets[i].r[2];
    }
    
    std::ofstream os;
    os.open("target.txt");
    for(int i = 0; i < numP; i++){
        os << i << "\t";
        os << targets[i].r[0] << "\t";
        os << targets[i].r[1] << "\t";
        os << targets[i].r[2] << std::endl;
    }  
}

void Model::calForcesHelper(int i, int j, double F[3]) {
    double r[dimP], dist;

    dist = 0.0;
    for (int k = 0; k < dimP; k++) {
        F[k] = 0.0;
        r[k] = (particles[j]->r[k] - particles[i]->r[k]) / radius;
        dist += pow(r[k], 2.0);
    }
    dist = sqrt(dist);
    if (dist < 2.0) {
        std::cout << "overlap " << i << "\t" << j << "\t"<< this->timeCounter <<std::endl;
        dist = 2.001;
    }
    if (dist < cutoff) {
        double Fpp = LJ * (-12.0 * pow((rm / dist), 12) / dist + 12.0 * pow((rm / dist), 7) / dist);
        Fpp += -9e-13 * exp(-33.3 * (dist - 2.0));
        for (int k = 0; k < dimP; k++) {
            F[k] = Fpp * r[k] / dist;

        }
    }


}

void Model::calForces() {
    double r[dimP], dist, F[3];
    for (int i = 0; i < numP; i++) {
        for (int k = 0; k < dimP; k++) {
            particles[i]->F[k] = 0.0;
        }
    }

    if(!cellListFlag){
    for (int i = 0; i < numP - 1; i++) {
        for (int j = i + 1; j < numP; j++) {
            calForcesHelper(i, j, F);
            for (int k = 0; k < dimP; k++) {
                particles[i]->F[k] += F[k];
                particles[j]->F[k] += -F[k];
            }
        }
    }
    } else{
        
        for (int i = 0; i < numP; i++) {
            std::vector<int> mapTable;
            mapTable.assign(numP,0);
            std::vector<int> nblist = 
            cellList->getNeighbors(particles[i]->r[0],particles[i]->r[1],particles[i]->r[2]);
            for (int j = 0; j < nblist.size(); j++){
                mapTable[nblist[j]] = 1;
                if (i!=nblist[j]){
                    calForcesHelper(i, nblist[j], F);
                    for (int k = 0; k < dimP; k++) {
                        particles[i]->F[k] += F[k];
                    }
                }
            }
/*           
            for (int j = 0; j < numP; j++){
                if (mapTable[j]!=1){
                        double r[3], dist;
                        dist = 0.0;
                        for (int k = 0; k < dimP; k++) {

                            r[k] = (particles[j]->r[k] - particles[i]->r[k]) / radius;
                            dist += pow(r[k], 2.0);
                        }
                        dist = sqrt(dist);
                        if (dist < cutoff) {
                            std::cout << "particle: " << j <<"uncovered! " << dist <<std::endl;
                            int idx1[3],idx2[3];
                            cellList->getParticleIdx(particles[i]->r[0],particles[i]->r[1],particles[i]->r[2],idx1);
                            cellList->getParticleIdx(particles[j]->r[0],particles[j]->r[1],particles[j]->r[2],idx2);
                            
                            std::cout << "particle i idx: " << idx1[0] << "\t" << idx1[1] << "\t"<<idx1[2]<<std::endl;
                            std::cout << "particle j idx: " << idx2[0] << "\t" << idx2[1] << "\t" <<idx2[2]<<std::endl;
                            cellList->printCellContent(idx1);
                            cellList->printCellContent(idx2);
                            cellList->buildList(particles);
                            cellList->printCellContent(idx1);
                            cellList->printCellContent(idx2);
                            
                        }   
                    
                }
            }         
            
            
            double F2[3], diff[3];
            diff[0]=0.0;
            diff[1]=0.0;
            diff[2]=0.0;
            for (int j = 0; j < numP; j++){
                if (i!=j){
                    calForcesHelper(i, j, F2);
                    for (int k = 0; k < dimP; k++) {
                        diff[k] += F2[k];
                    }
                    
                }
            }              
            
            std::cout << "particle: " << i << std::endl; 
                    std::cout << diff[0] - particles[i]->F[0] << std::endl;                    
                    std::cout << diff[1] - particles[i]->F[1] << std::endl;
                    std::cout << diff[2] - particles[i]->F[2] << std::endl;
            }
*/             
        }    
    }
    
}
    


void Model::createInitialState(){

    this->readxyz(iniFile);
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
        for (int j = 0; j < dimP; j++){
            os << particles[i]->r[j]/radius << "\t";
        }
        os << particles[i]->phi<< "\t";
        os << particles[i]->theta<< "\t";
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
        linestream >> particles[i]->r[2];
        linestream >> particles[i]->phi;
        linestream >> particles[i]->theta;
    }
    for (int i = 0; i < numP; i++) {
        particles[i]->r[0] *=radius;
        particles[i]->r[1] *=radius;
        particles[i]->r[2] *=radius;
        particles[i]->ori_vec[0][0] = cos(particles[i]->phi)*sin(particles[i]->theta);
        particles[i]->ori_vec[0][1] = sin(particles[i]->phi)*sin(particles[i]->theta);
        particles[i]->ori_vec[0][2] = cos(particles[i]->theta);
    }
    this->updateBodyFrameVec();
    
    if (cellListFlag){
        cellList->buildList(particles);
    }
    
    is.close();
}

void Model::updateBodyFrameVec(){
    
    double thresh = 0.99999;
    double norm;
    
      for (int ii = 0; ii < numP; ii++){
          // set n2, n3 to zero
          for (int i = 0; i < 3; i++){
              particles[ii]->ori_vec[1][i] = 0.0;
              particles[ii]->ori_vec[2][i] = 0.0;
              
          }
          
          norm = 0.0;
          for(int i = 0; i < 3; i++){
              norm = norm + pow(particles[ii]->ori_vec[0][i],2.0);
          }
                    norm = sqrt(norm);
          for(int i = 0; i < 3; i++){
              particles[ii]->ori_vec[0][i] /= norm;
          }          
          
          
          // first do some safe-guard to prevent numerical instability
          if(particles[ii]->ori_vec[0][2] >= thresh){
              particles[ii]->ori_vec[0][0] = 0.0;
              particles[ii]->ori_vec[0][1] = 0.0;
              particles[ii]->ori_vec[0][2] = 1.0;
          }
      
          if(particles[ii]->ori_vec[0][2] <= -thresh){
              particles[ii]->ori_vec[0][0] = 0.0;
              particles[ii]->ori_vec[0][1] = 0.0;
              particles[ii]->ori_vec[0][2] = -1.0;
          }
      

           
          // first consider degenrate case that n1 = e_z
          if(abs(particles[ii]->ori_vec[0][2]) >= thresh){
             particles[ii]->ori_vec[1][0]=0;
             particles[ii]->ori_vec[1][1]=1.0;
             particles[ii]->ori_vec[1][2]=0;
          }
          else {
              for (int i = 0; i < 3; i++){
//                  for (int j = 0; j < 3; j++){
                      for (int k = 0; k < 3; k++){
                          particles[ii]->ori_vec[1][i] += eps[i][2][k]*particles[ii]->ori_vec[0][k];
                      
                      
                      }
//                  }
            }
            norm = 0.0;
            for (int i = 0; i < 3; i++) {
                norm = norm + pow(particles[ii]->ori_vec[1][i],2.0);
            }
            norm = sqrt(norm);
            for (int i = 0; i < 3; i++) {
                particles[ii]->ori_vec[1][i] /= norm;
            }
          }
          for (int i = 0; i < 3; i++){
                  for (int j = 0; j < 3; j++){
                      for (int k = 0; k < 3; k++){
                          particles[ii]->ori_vec[2][i] += 
                                  eps[i][j][k]*particles[ii]->ori_vec[1][j]
                                  *particles[ii]->ori_vec[0][k];
                      
                      
                      }
                  }
             }
                
    //                for(int i = 0; i < 3; i++){
    //                    for(int j = 0; j <3;j++){
    //                        std::cout << particles[ii]->ori_vec[i][j] << std::endl;
                        
    //                    }
    //                }
     }
}

void Model::getPermutator() {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                eps[i][j][k] == 0;
                if (i == 0 && j == 1 && k == 2) {
                    eps[i][j][k] = 1;
                }
                if (i == 1 && j == 2 && k == 0) {
                    eps[i][j][k] = 1;
                }
                if (i == 2 && j == 0 && k == 1) {
                    eps[i][j][k] = 1;
                }
                if (i == 2 && j == 1 && k == 0) {
                    eps[i][j][k] = -1;
                }
                if (i == 0 && j == 2 && k == 1) {
                    eps[i][j][k] = -1;
                }
                if (i == 1 && j == 0 && k == 2) {
                    eps[i][j][k] = -1;
                } 
            }
        }
    }
}