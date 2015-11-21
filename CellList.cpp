#include "CellList.h"

CellList::CellList(double cutDist, int dim0, int maxCount0, int box_x0,
        double del_x0, int box_y0, double del_y0, int box_z0, double del_z0) :
cutDistance(cutDist), dim(dim0), maxCount(maxCount0), boxSize_x(box_x0),
boxSize_y(box_y0), boxSize_z(box_z0), del_x(del_x0), del_y(del_y0), del_z(del_z0) {

    nbin_x = (int) (boxSize_x / del_x) + 1;
    nbin_y = (int) (boxSize_y / del_y) + 1;
    nbin_z = (int) (boxSize_z / del_z) + 1;
    if (dim == 2) {
        nbin_z = 1;
    }
    nbin = nbin_x * nbin_y * nbin_z;
    min_x = -0.5 * boxSize_x;
    min_y = -0.5 * boxSize_y;
    min_z = -0.5 * boxSize_z;
    this->setup();
}

CellList::Idx_3d CellList::coordToIdx(double x, double y, double z) {

    int z_bin, x_bin, y_bin;

    if (dim == 2) {
        z_bin = 0;
        x_bin = (int) ((x - min_x) / del_x);
        y_bin = (int) ((y - min_y) / del_y);
    } else {
        x_bin = (int) ((x - min_x) / del_x);
        y_bin = (int) ((y - min_y) / del_y);
        z_bin = (int) ((z - min_z) / del_z);
    }
    Idx_3d idx(z_bin, y_bin, x_bin);
    return idx;
}

std::vector<int> CellList::getNeighbors(double x, double y, double z) {

    Idx_3d idx = coordToIdx(x,y,z);
    std::vector<int> res;
    int idx_x = std::get<2>(idx);
    int idx_y = std::get<1>(idx);
    int idx_z = std::get<0>(idx);
    int nb_idx = (*oneDIdx)[idx_z][idx_y][idx_x];
    for (auto& cellIdx : cellNeighborIdxList[nb_idx]) {
        int idx_x = std::get<2>(cellIdx);
        int idx_y = std::get<1>(cellIdx);
        int idx_z = std::get<0>(cellIdx);
        int count = (*cellListCount)[idx_z][idx_y][idx_x];
        for (int i = 0; i < count; i++) {
            res.push_back((*cellList)[idx_z][idx_y][idx_x][i]);
        }
    }
    return res;
}
// build cellList

void CellList::buildList(const Model::state &s) {
    int totalCount = 0;
    for (int i = 0; i < s.size(); i++) {
        int count = 0; 
        Idx_3d idx = coordToIdx(s[i]->r[0], s[i]->r[1], s[i]->r[2]);
        threeDIdx.push_back(idx);
        int idx_x = std::get<2>(idx);
        int idx_y = std::get<1>(idx);
        int idx_z = std::get<0>(idx);
        (*cellList)[idx_z][idx_y][idx_x][count] = i;        
        count++;
        (*cellListCount)[idx_z][idx_y][idx_x]=count;
        totalCount++;
    }
    std::cout <<totalCount<< "particle binned" << std::endl;
}

void CellList::setup() {

    std::array<Array4D_type::index, 4> dim1 = {nbin_z, nbin_y, nbin_x, maxCount};
    std::array<Array4D_type::index, 3> dim2 = {nbin_z, nbin_y, nbin_x};
    cellList = std::make_shared<Array4D_type>(dim1);
    oneDIdx = std::make_shared<Array3D_type>(dim2);
    cellListCount = std::make_shared<Array3D_type>(dim2);

    int count = 0;
    // now build the cellNeighborIdxList
    for (int i = 1; i < nbin_x - 1; i++) {
        for (int j = 1; j < nbin_y - 1; j++) {
            if (dim == 2) {
                cellNeighborIdxList.push_back(std::list<Idx_3d>());
                for (int ii = -1; ii < 2; ii++) {
                    for (int jj = -1; jj < 2; jj++) {
                        int idx_x = i + ii;
                        int idx_y = j + jj;
                        int idx_z = 0;
                        Idx_3d idx(idx_z, idx_y, idx_x);
                        (*oneDIdx)[0][j][i] = count;
                        cellNeighborIdxList[count].push_back(idx);
                        
                    }
                }
                count++;
            } else {

                for (int k = 1; k < nbin_z - 1; k++) {
                    cellNeighborIdxList.push_back(std::list<Idx_3d>());
                    for (int ii = -1; ii < 2; ii++) {
                        for (int jj = -1; jj < 2; jj++) {
                            for (int kk = -1; kk < 2; kk++) {
                                int idx_x = i + ii;
                                int idx_y = j + jj;
                                int idx_z = k + kk;
                                Idx_3d idx(idx_z, idx_y, idx_x);
                                (*oneDIdx)[k][j][i] = count;
                                cellNeighborIdxList[count].push_back(idx);

                            }
                        }
                    }
                }
                count++;
            }
        }
    }

}

void CellList::printCellList() const{
    for(int nb_idx = 0; nb_idx < nbin; nb_idx++){
        std::cout << nb_idx << std::endl;
        for (auto& cellIdx : cellNeighborIdxList[nb_idx]) {
            int idx_x = std::get<2>(cellIdx);
            int idx_y = std::get<1>(cellIdx);
            int idx_z = std::get<0>(cellIdx);
            int count = (*cellListCount)[idx_z][idx_y][idx_x];
            std::cout << idx_x << "\t" << idx_y << "\t" << idx_z << std::endl;
        }
    }
}


void CellList::printParticleList() const{
    for(int i = 0; i < threeDIdx.size(); i++){
        Idx_3d cellIdx = threeDIdx[i];
        int idx_x = std::get<2>(cellIdx);
        int idx_y = std::get<1>(cellIdx);
        int idx_z = std::get<0>(cellIdx);
        std::cout << i<< "\t" <<idx_x << "\t" << idx_y << "\t" << idx_z << std::endl;
       
    }
}