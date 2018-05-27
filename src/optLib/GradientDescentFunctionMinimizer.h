#pragma once

#include "ObjectiveFunction.h"

#include <string>
#include <iostream>
#include <cmath>
#include <cfloat>
//#define DEBUG

class GradientDescentFunctionMinimizer
{
  public:
	GradientDescentFunctionMinimizer(int maxIterations = 100, double solveResidual = 1e-5, int maxLineSearchIterations = 15)
		: maxIterations(maxIterations), solveResidual(solveResidual), maxLineSearchIterations(maxLineSearchIterations)
	{
	}

	virtual ~GradientDescentFunctionMinimizer() {}

	int getLastIterations() { return lastIterations; }

	virtual bool minimize(ObjectiveFunction *function, VectorXd &x)
	{

		//number of parameters...
		int N = (int)x.size();
		resize(xi, N);
		resize(dx, N);
		resize(gradient, N);

		xi = x;

		bool optimizationConverged = false;

		int i = 0;
		for (; i < maxIterations; i++)
		{
			computeSearchDirection(function, xi, dx);
			double actualPerf = dx.norm();
#ifdef DEBUG_OLD
			std::cout << "N°Iteration (minimization) " << i << " actualPerf (minimization) " << actualPerf << std::endl;
#endif

			if (actualPerf < solveResidual)
			{
				optimizationConverged = true;
				break;
			}

			doLineSearch(function, dx, xi);
		}

		lastIterations = i;

		//p now holds the parameter values at the start of the iteration...
		x = xi;

		//and done!
		return optimizationConverged;
	}

  protected:
	// Since the gradient of a function gives the direction of steepest descent, all one needs to do is go in that direction...
	virtual void computeSearchDirection(ObjectiveFunction *function, const VectorXd &x, VectorXd &dx)
	{
#ifdef DEBUG
		std::cout << "\n> DEBUG computeSearchDirection - size of x (input) : " << x.rows() << " by " << x.cols() << std::endl;
		std::cout << "> DEBUG computeSearchDirection - size of dx (input) : " << dx.rows() << " by " << dx.cols() << std::endl;
#endif
		// Ex. 1.1
		//VectorXd grad;
		//grad.resize(x.rows(), 1);
		dx.setZero();
		function->addGradientTo(dx, x);

		//Store it
		//dx = - grad;
		dx *= -1;
		//dx[0] = -grad[0];
		//dx[1] = -grad[1];
		
		//Note : no normalization, otherwise we loose information ?
	}

	virtual void doLineSearch(ObjectiveFunction *function, const VectorXd &dx, VectorXd &xi)
	{
#ifdef DEBUG
		std::cout << "\n> DEBUG doLineSearch - size of x (input) : " << x.rows() << " by " << x.cols() << std::endl;
		std::cout << "> DEBUG doLineSearch - size of dx (input) : " << dx.rows() << " by " << dx.cols() << std::endl;
#endif
		// Ex. 1.1
		bool isInferior = false;
		VectorXd tmpx;
		tmpx.resize(xi.rows(), 1);
		double refValue = function->computeValue(xi);
		double currentAlpha = alpha;

		tmpx = xi + dx*currentAlpha;
		double tmpValue = function->computeValue(tmpx);

		for (int i = 0; i < maxLineSearchIterations || !std::isfinite(tmpValue); i++) //Can not finnish
		{

			if (tmpValue < refValue && std::isfinite(tmpValue))
			{ // We have a lower point, we stop here
				break;
			}

			//Next iteration
			currentAlpha *= beta;
			tmpx = xi + dx*currentAlpha;
			tmpValue = function->computeValue(tmpx);
		}

		//Storage
		xi = tmpx;
	}

  protected:
	double solveResidual = 1e-5;
	const int maxIterations = 100;
	const int maxLineSearchIterations = 15;

	const double alpha = 1;
	const double beta = 0.5;

	VectorXd xi, dx, gradient;

	// some stats about the last time `minimize` was called
	int lastIterations = -1;
};
