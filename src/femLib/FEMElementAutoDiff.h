#pragma once

#include "FEMElement.h"

#include <AutoDiff.h>

template <class S> using Vector2 = Eigen::Matrix<S, 2, 1>;
template <class S> using Matrix2 = Eigen::Matrix<S, 2, 2>;

typedef AutoDiff<double, double> AD;
typedef AutoDiff<AD, AD> ADD;

/* This is the code from the tutorial using this AutoDiff implementation:
 *
 *		// df/dx1
		AD x1 = M_PI; x1.deriv() = 1;
		AD x2 = 3; x2.deriv() = 0;
		AD f = sin(x1*x2);
		double f_val = f.value();
		double dfdx1 = f.deriv();

		// ddf/dx1dx2
		ADD x1 = M_PI;
		x1.value().deriv() = 1; x1.deriv().value() = 0;
		ADD x2 = 3;
		x2.value().deriv() = 0; x1.deriv().value() = 1;
		ADD f = sin(x1*x2);
		double f_val = f.value().value();
		double ddfdx1dx2 = f.deriv().deriv();

		// ddf/dx2dx2
		ADD x1 = M_PI;
		x1.value().deriv() = 0; x1.deriv().value() = 0;
		ADD x2 = 3;
		x2.value().deriv() = 1; x1.deriv().value() = 1;
		ADD f = sin(x1*x2);
		double f_val = f.value().value();
		double ddfdx1dx2 = f.deriv().deriv();
*/

/**
	This class implements Constant Strain Triangles elements in 2D
*/
class FEMElementAutoDiff : public FEMElement {
public:
	FEMElementAutoDiff(const std::array<int, 3> &nodeIndices, const VectorXd &X)
		: FEMElement(nodeIndices, X) {
	}
	~FEMElementAutoDiff(){
	}

	virtual double getEnergy(const VectorXd& x, const VectorXd& X){

		// Ex 2.3

	}

	virtual void addEnergyGradientTo(const VectorXd& x, const VectorXd& X, VectorXd& grad) {

		// Ex. 2.3

	}

	virtual void addEnergyHessianTo(const VectorXd& x, const VectorXd& X, std::vector<Tripletd>& hesEntries) {

		// Ex. 2.3

	}

private:
	// This is the same function as FEMElement::computeGradient,
	// however with a template parameter T.
	template<class T>
	void computeDeformationGradientT(const Vector2<T> (&x)[3], Matrix2<T>& dxdX) const {
		//edge vectors
		Vector2<T> v1 = x[1] - x[0];
		Vector2<T> v2 = x[2] - x[0];
		Matrix2<T> dx;
		dx << v1[0], v2[0],
				v1[1], v2[1];
		Matrix2<T> dxInvT;
		for (int i = 0; i < 4; ++i)
			dxInvT.data()[i] = dXInv.data()[i];
		dxdX = dx * dxInvT;
	}
};
