#pragma once

#include "ObjectiveFunction.h"
#include "GradientDescentFunctionMinimizer.h"

//#define DEBUG

/**
	use Newton's method to optimize a function. p will store the final value that minimizes the function, and its initial value
	is used to start the optimization method.

	Task: find p that minimize f(p). This means that df/dp(p) = 0.
	df/dp(p+dp) ~ df/dp(p) + d/dp(df/dp) * dp = 0 ==> -df/dp(p) = d/dp(df/dp) * dp
	Iterating the above, will hopefully get p that minimizes f.
*/
class NewtonFunctionMinimizer : public GradientDescentFunctionMinimizer
{
  public:
	NewtonFunctionMinimizer(int maxIterations = 100, double solveResidual = 0.0001, int maxLineSearchIterations = 15)
		: GradientDescentFunctionMinimizer(maxIterations, solveResidual, maxLineSearchIterations) {}

	virtual ~NewtonFunctionMinimizer() {}

  protected:
	// The search direction is given by -Hinv * g, which is not optimal, so we solve a system
	virtual void computeSearchDirection(ObjectiveFunction *function, const VectorXd &x, VectorXd &dx)
	{
#ifdef DEBUG
		std::cout << "\n> DEBUG computeSearchDirection - size of x (input) : " << x.rows() << " by " << x.cols() << std::endl;
		std::cout << "\n> DEBUG computeSearchDirection - size of dx (input) : " << dx.rows() << " by " << dx.cols() << std::endl;
#endif
		// Ex 1.3

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		Eigen::SparseMatrix<double> hessianSparse(x.rows(), x.rows());
		Eigen::SparseMatrix<double> rightSide(x.rows(), 1);
		Eigen::SparseMatrix<double> solved(x.rows(), 1);
		std::vector<Tripletd> hessianEntries;
		VectorXd grad;

		grad.resize(x.rows(), 1);
		grad.setZero();

		function->addHessianEntriesTo(hessianEntries, x);
		hessianSparse.setFromTriplets(hessianEntries.begin(), hessianEntries.end());
		//hessianSparse *= 0.0001;

		function->addGradientTo(grad, x);
		grad = -1 * grad;
		rightSide = grad.sparseView();

		solver.compute(hessianSparse);
		solved = solver.solve(rightSide);

		dx = MatrixXd(solved);

/*
total energy = 122.621
# iterations = 100000
min def at   = 3.84296e-10
max def at   = 0.0428256
*/

#ifdef DEBUG
		std::cout << "\n> DEBUG computeSearchDirection - hessianSparse (inside) : " << hessianSparse.rows() << " by " << hessianSparse.cols() << std::endl;
		MatrixXd printableHessian = MatrixXd(hessianSparse);
		std::cout << printableHessian << std::endl;
		std::cout << "\n> DEBUG computeSearchDirection - rightSide (input) : " << rightSide.rows() << " by " << rightSide.cols() << std::endl;
		MatrixXd printableGrad = MatrixXd(rightSide);
		std::cout << printableGrad << std::endl;
#endif
	}

  public:
	SparseMatrixd H;
	std::vector<Tripletd> hessianEntries;
};
