#pragma once
#include <vector>
#include <list>
#include <tuple>
#include "boost/multi_array.hpp"
#include "model.h"


class CellList{
public:
	CellList(double cutDist, int dim0, int maxCount0, int box_x,
        double del_x, int box_y, double del_y, int box_z, double del_z);
	~CellList(){}
	std::vector<int> getNeighbors(double x, double y, double z);
	void buildList(const Model::state &s);
	void setup();
        void printCellList() const;
        void printParticleList() const;
private:
    typedef std::tuple<int,int,int> Idx_3d;
	int dim, maxCount;
	double boxSize_x, boxSize_y, boxSize_z;
        double min_x, min_y, min_z;
	double del_x, del_y, del_z;
	int nbin_x, nbin_y, nbin_z, nbin;
	double cutDistance;
	Idx_3d coordToIdx(double x, double y, double z);

	typedef boost::multi_array<int, 4> Array4D_type;
        typedef boost::multi_array<int, 3> Array3D_type;
	std::shared_ptr<Array4D_type> cellList;
	std::shared_ptr<Array3D_type> cellListCount, oneDIdx ;
//	boost::multi_array<int,3> cellListCount;
//	std::vector<std::vector<int>> cellList;
	std::vector<std::list<Idx_3d>> cellNeighborIdxList;
        std::vector<Idx_3d> threeDIdx;
};