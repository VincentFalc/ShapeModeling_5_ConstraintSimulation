#pragma once

#include "Element.h"
#include "TotalEnergyFunction.h"

/**
	This class implements a generic sim mesh for deformable objects: collection of nodes connected to each other using different types of elements
*/
class SimulationMesh {

public:
	SimulationMesh();
	~SimulationMesh();

	// create and load meshes
	struct TriangleMesh {
		std::vector<Vector2d> v;				// vertices
		std::vector<std::array<int, 3>> tri;	// triangles
	};
	static void generateSquareTriMesh(TriangleMesh &triMesh, double startX=-4.5, double startY=0, double dX=1, double dY=1, int xSize=10, int ySize=10);
	void loadTriangleMesh(const TriangleMesh &triMesh);
	void createSpringsFromTriangles();
	template<class ElementType>	void createElementsFromTriangles();

	// add external forces
	void addGravityForces(const Vector2d &g);
	void applyForceAt(int i, const Vector2d &force);

	// solve
	template<class Minimizer> void solve_statics();
	void solveGradientDescent();
	void solveNewtonsMethod();

	// used for the report
	void testGradient();
	void testHessian();

	// pinned nodes
	void togglePinnedNode(int i);
	void setPinnedNode(int nodeIndex, const Vector2d& point);
	void unpinNode(int nodeIndex);

	void computeDefoEnergyPerNode(Eigen::VectorXd &defEnergyPerNode);

	void clear();

public:
	// number of nodes
	int numNodes, numElements;

	// for each node we will store the position and velocity, rest configuration, mass and external forces acting on it
	VectorXd x, m, f_ext, X, xSolver;

	// list of elements that connects the nodes to each other
	std::vector<std::array<int, 3>> triangles;

	// a list of temporary used to pin points to locations that are fixed in space...
	std::vector<Element*> pinnedNodeElements;

	//this is the objective function that we use for simulations...
	TotalEnergyFunction* energyFunction;
};



template<class ElementType>
void SimulationMesh::createElementsFromTriangles()
{
	energyFunction->elements.clear();
	for (int i=0; i<numElements; i++){
		Element* newElem = new ElementType(triangles[i], X);
		energyFunction->elements.push_back(newElem);
	}

	// update mass
	m.setZero();
	for(auto e : energyFunction->elements) {
		for (int i = 0; i < e->getNumNodes(); ++i) {
			for (int j=0; j<2; j++)
				this->m[2*e->getNodeIndex(i) + j] += e->getMass() / (double)e->getNumNodes();
		}
	}
}

template<class Minimizer>
void SimulationMesh::solve_statics(){
	xSolver = x;

	double functionValue = energyFunction->computeValue(xSolver);

	Minimizer minimizer(50);
	minimizer.minimize(energyFunction, xSolver);

	x = xSolver;
}
