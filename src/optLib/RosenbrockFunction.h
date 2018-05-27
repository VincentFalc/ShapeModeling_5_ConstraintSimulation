#pragma once

#include "ObjectiveFunction.h"
//#define DEBUG

class RosenbrockFunction : public ObjectiveFunction
{
  public:
	RosenbrockFunction()
	{
		a = 1;
		b = 100;
	}

	virtual double computeValue(const VectorXd &x)
	{

		// Ex 1.1
		// return f(x)
		double xV = x[0];
		double yV = x[1];

#ifdef DEBUG
		std::cout << "\n> DEBUG computeValue - size of x (input) : " << x.rows() << " by " << x.cols() << std::endl;
		std::cout << "> DEBUG computeValue - values of x (input) : " << xV << " & " << yV << std::endl;
		std::cout << "> DEBUG computeValue - value of f(x) (output) : " << (a - xV) * (a - xV) + b * (yV - (xV * xV)) * (yV - (xV * xV)) << std::endl;
#endif

		return (a - xV) * (a - xV) + b * (yV - (xV * xV)) * (yV - (xV * xV));
	}

	virtual void addGradientTo(VectorXd &grad, const VectorXd &x)
	{
#ifdef DEBUG
		std::cout << "\n> DEBUG addGradientTo - size of x (input) : " << x.rows() << " by " << x.cols() << std::endl;
		std::cout << "> DEBUG addGradientTo - size of grad (input) : " << grad.rows() << " by " << grad.cols() << std::endl;
#endif
		assert(x.rows()==2);
		assert(grad.rows()==2);

		double xV = x[0];
		double yV = x[1];

		//Calculate the directionnal derivatives
		double dfdx = 2 * (-a + xV + 2 * b * xV * (xV*xV - yV));
		double dfdy = 2 * b * (yV - xV*xV);

		assert(std::isfinite(dfdx));
		assert(std::isfinite(dfdy));

		//Store it
		grad[0] += dfdx; //+ grad[0];
		grad[1] += dfdy; //+ grad[1]; /// ? Adds ? Or replace ? Or appends ?
	}

	virtual void addHessianEntriesTo(std::vector<Tripletd> &hessianEntries, const VectorXd &x)
	{
#ifdef DEBUG
		std::cout << "\n> DEBUG addHessianEntriesTo - size of x (input) : " << x.rows() << " by " << x.cols() << std::endl;
		std::cout << "> DEBUG addHessianEntriesTo - size of hessianEntries (input) : " << hessianEntries.size() << std::endl;
#endif
		assert(x.rows()==2);
		assert(hessianEntries.size()==0);

		// Ex 1.2
		// write d^2f/dx^2 in `hessianEntries`
		double xV = x[0];
		double yV = x[1];

		//Empty the Tripletd ?
		hessianEntries.push_back(Tripletd(0, 0, -4 * b * (yV - xV * xV) + 8 * b * xV * xV + 2));
		hessianEntries.push_back(Tripletd(0, 1, -4 * b * xV));
		hessianEntries.push_back(Tripletd(1, 0, -4 * b * xV));
		hessianEntries.push_back(Tripletd(1, 1, 2 * b));
#ifdef DEBUG
		std::cout << "> DEBUG addHessianEntriesTo - hessianEntries (inside) : " << hessianEntries.size() << std::endl;
		std::cout << "> DEBUG addHessianEntriesTo - hessianEntries (inside) : VAL1 = " << -4 * b * (yV - xV * xV) + 8 * b * xV * xV + 2 << std::endl;
		std::cout << "> DEBUG addHessianEntriesTo - hessianEntries (inside) : VAL2&3 =" << -4 * b * xV << std::endl;
		std::cout << "> DEBUG addHessianEntriesTo - hessianEntries (inside) : VAL4 =" << 2 * b << std::endl;
#endif
	}

	double a, b;
};
