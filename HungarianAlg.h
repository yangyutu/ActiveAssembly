#include <vector>
#include <iostream>
#include <limits>
#include <time.h>
// http://community.topcoder.com/tc?module=Static&d1=tutorials&d2=hungarianAlgorithm
using namespace std;
class AssignmentProblemSolver
{
private:
	// --------------------------------------------------------------------------
	// Computes the optimal assignment (minimum overall costs) using Munkres algorithm.
	// --------------------------------------------------------------------------
	void assignmentoptimal(int *assignment, long *cost, long *distMatrix, int nOfRows, int nOfColumns);
	void buildassignmentvector(int *assignment, bool *starMatrix, int nOfRows, int nOfColumns);
	void computeassignmentcost(int *assignment, long *cost, long *distMatrix, int nOfRows);
	void step2a(int *assignment, long *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
	void step2b(int *assignment, long *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
	void step3 (int *assignment, long *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
	void step4 (int *assignment, long *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim, int row, int col);
	void step5 (int *assignment, long *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
	// --------------------------------------------------------------------------
	// Computes a suboptimal solution. Good for cases with many forbidden assignments.
	// --------------------------------------------------------------------------
	void assignmentsuboptimal1(int *assignment, long *cost, long *distMatrixIn, int nOfRows, int nOfColumns);
	// --------------------------------------------------------------------------
	// Computes a suboptimal solution. Good for cases with many forbidden assignments.
	// --------------------------------------------------------------------------
	void assignmentsuboptimal2(int *assignment, long *cost, long *distMatrixIn, int nOfRows, int nOfColumns);
public:
	enum TMethod { optimal, many_forbidden_assignments, without_forbidden_assignments };
	AssignmentProblemSolver();
	~AssignmentProblemSolver();
	long Solve(vector<vector<long>>& DistMatrix,vector<int>& Assignment,TMethod Method=optimal);
};