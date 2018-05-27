#pragma once

#include "FEMElement.h"

/**
	This class implements Constant Strain Triangles elements in 2D
*/
class FEMElementAnalytic : public FEMElement {
public:
	FEMElementAnalytic(const std::array<int, 3> &nodeIndices, const VectorXd &X)
		: FEMElement(nodeIndices, X) {
	}
	~FEMElementAnalytic(){
	}

	virtual double getEnergy(const VectorXd& x, const VectorXd& X){
		// Ex. 2.3 (same as FEMElementFD
		return 0;
	}

	virtual void addEnergyGradientTo(const VectorXd& x, const VectorXd& X, VectorXd& grad) {
		// Ex. 2.3
	}

	virtual void addEnergyHessianTo(const VectorXd& x, const VectorXd& X, std::vector<Tripletd>& hesEntries) {
		// Ex. 2.3
	}
};
