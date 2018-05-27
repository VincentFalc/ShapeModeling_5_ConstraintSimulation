#pragma once

#include "FEMElement.h"
//#define DEBUG
//#define DEBUG_Verbose

/**
	This class implements Constant Strain Triangles elements in 2D
*/
class FEMElementFD : public FEMElement
{
  public:
	FEMElementFD(const std::array<int, 3> &nodeIndices, const VectorXd &X)
		: FEMElement(nodeIndices, X)
	{
	}
	~FEMElementFD()
	{
	}

	virtual double getEnergy(const VectorXd &x, const VectorXd &X)
	{
		// Ex 2.1
		Vector2d xTMP[3];
		xTMP[0] = getNodePos(0, x);
		xTMP[1] = getNodePos(1, x);
		xTMP[2] = getNodePos(2, x);
		Matrix2d dxdX; //2*2 ?

		//computeDeformationGradient(const Vector2d (&x)[3], Matrix2d& dxdX)
		computeDeformationGradient(xTMP, dxdX);

#ifdef DEBUG_Verbose
		std::cout << "Energy calculus - dxdX: \t" << dxdX << std::endl;
#endif

		MatrixXd C = (dxdX.transpose() * dxdX);
		double logdetF = log(dxdX.determinant());

#ifdef DEBUG_Verbose
		std::cout << "Energy calculus - dxdX.determinant(): \t" << dxdX.determinant() << std::endl;
		std::cout << "Energy calculus - logdetF: \t" << logdetF << std::endl;
#endif

		double Energy = (shearModulus * 0.5 * (C.trace() - 2) - shearModulus * logdetF + bulkModulus * 0.5 * logdetF * logdetF);

#ifdef DEBUG
		std::cout << "Energy calculus : \t" << Energy << std::endl;
		std::cout << "==== End of energy calculation ====" << std::endl;
#endif

		/*
		if(!std::isfinite(Energy)){
			std::cout << "There is a problem" << std::endl;
		}
*/

		return restShapeArea * Energy;
	}

	virtual void addEnergyGradientTo(const VectorXd &x, const VectorXd &X, VectorXd &grad)
	{
		// Ex 2.1

		//For each point, we move it a bit in a +eps and a -eps, on each coordinate, to compute the finite difference
		int indexPoint[3];
		indexPoint[0] = getNodeIndex(0);
		indexPoint[1] = getNodeIndex(1);
		indexPoint[2] = getNodeIndex(2);

		VectorXd modified_x_A;
		VectorXd modified_x_B;

		for (int j = 0; j < 2; j++)
		{
			for (int i = 0; i < 3; i++)
			{
				modified_x_A = x;
				modified_x_B = x;

				modified_x_A[2 * indexPoint[i] + j] += h;
				modified_x_B[2 * indexPoint[i] + j] -= h;
				double Etmp = (getEnergy(modified_x_A, X) - getEnergy(modified_x_B, X)) / (2 * h);
#ifdef DEBUG
				std::cout << "Energy : " << Etmp << std::endl;
#endif
				grad[2 * indexPoint[i] + j] += Etmp;
			}
		}
}

virtual void addEnergyHessianTo(const VectorXd &x, const VectorXd &X, std::vector<Tripletd> &hesEntries)
{
	// Ex 2.2

	//Liste of necessary energies
	VectorXd modified_x_A_p_B_p;
	VectorXd modified_x_A_p_B_m;
	VectorXd modified_x_A_m_B_p;
	VectorXd modified_x_A_m_B_m;

	double Energy_A_p_B_p;
	double Energy_A_p_B_m;
	double Energy_A_m_B_p;
	double Energy_A_m_B_m;

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			//We get the indices of two nodes
			int nodeAIndex = getNodeIndex(i);
			int nodeBIndex = getNodeIndex(j);

			//For the corner of this hessian bloc : xx, xy, yx, yy
			for (int x_val = 0; x_val < 2; x_val++) //Refers to the verical position in the actual 4-sized-square in the Hessian matrix
			{
				for (int y_val = 0; y_val < 2; y_val++) //Refers to the horizontal position in the actual 4-sized-square in the Hessian matrix
				{
					//Set the intiales values of the
					modified_x_A_p_B_p = x;
					modified_x_A_p_B_m = x;
					modified_x_A_m_B_p = x;
					modified_x_A_m_B_m = x;

					//Modify positions
					modified_x_A_p_B_p[2 * nodeAIndex + x_val] += h; //Ax or Ay
					modified_x_A_p_B_p[2 * nodeBIndex + y_val] += h; //Bx or By

					modified_x_A_p_B_m[2 * nodeAIndex + x_val] += h;
					modified_x_A_p_B_m[2 * nodeBIndex + y_val] -= h;

					modified_x_A_m_B_p[2 * nodeAIndex + x_val] -= h;
					modified_x_A_m_B_p[2 * nodeBIndex + y_val] += h;

					modified_x_A_m_B_m[2 * nodeAIndex + x_val] -= h;
					modified_x_A_m_B_m[2 * nodeBIndex + y_val] -= h;

					//Energy calculations
					Energy_A_p_B_p = getEnergy(modified_x_A_p_B_p, X);
					Energy_A_p_B_m = getEnergy(modified_x_A_p_B_m, X);
					Energy_A_m_B_p = getEnergy(modified_x_A_m_B_p, X);
					Energy_A_m_B_m = getEnergy(modified_x_A_m_B_m, X);

					//We compute the Hessian values (Gradients inside)
					double Hess_A_B_x_y = (((Energy_A_p_B_p - Energy_A_p_B_m) / (2 * h)) - ((Energy_A_m_B_p - Energy_A_m_B_m) / (2 * h))) / (2 * h);

					//We add the value in the hessian
					hesEntries.push_back(Tripletd(2 * nodeAIndex + x_val, 2 * nodeBIndex + y_val, Hess_A_B_x_y));
				}
			}

		}
	}

}

private:
double h = 1e-8; // step size
};
