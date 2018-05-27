#pragma once
//#define DEBUG

#include "Element.h"

/**
	This class implements the interface for an elementary energy unit. As a function of deformed, undeformed,
	and other parameters, such as boundary conditions, each class that extends this one will define a potential energy.
	The deformed energy depends on a number of nodes.
*/
class Spring : public Element {

public:
	Spring(const std::array<int, 2> &nodeIndices, const VectorXd &X)
		: nodeIndices(nodeIndices) {
	}
	virtual ~Spring() {}

	// Returns the number of nodes this unit depends on
	virtual int getNumNodes() const {
		return 2;
	}
	// Returns the global index of node `i`
	virtual int getNodeIndex(int i) const {
		return nodeIndices[i];
	}

	// Returns the element's mass
	virtual double getMass() const {
		return 0;
	}

	virtual double getLength(const Vector2d & Pos_0, const Vector2d & Pos_1){
		double delta_x = Pos_0[0]-Pos_1[0];
		double delta_y = Pos_0[1]-Pos_1[1];
		//double l = sqrt(delta_x*delta_x + delta_y*delta_y); // sqrt( (x1-x2)² + (y1-y2)² )
		return sqrt(delta_x*delta_x + delta_y*delta_y);  // sqrt( (x1-x2)² + (y1-y2)² )
	}

	virtual double getNorm(const Vector2d & Pos_Delta){
		//double l = sqrt(delta_x*delta_x + delta_y*delta_y); // sqrt( (x1-x2)² + (y1-y2)² )
		return sqrt(Pos_Delta[0]*Pos_Delta[0] + Pos_Delta[1]*Pos_Delta[1]);  // sqrt( (x1-x2)² + (y1-y2)² )
	}

	virtual Vector2d computeGrad(const double k, const Vector2d l, const double  L_norm, Vector2d& u){
		double l_norm = getNorm(l);
		double cste = k * ((l_norm-L_norm)/L_norm);
		assert(std::isfinite(cste));

		Vector2d grad(2);
		grad.setZero();

		grad[0] += cste * u[0];
		grad[1] += cste * u[1];

		return grad;
	}

	virtual double computeEnergy(double elongation){
		return 0.5*k*elongation*elongation;
	}

	virtual double computeEnergy(Vector2d l, double L_length){
		double l_length = getNorm(l);
		return 0.5 * k * ((l_length-L_length)/L_length) * ((l_length-L_length)/L_length);
	}

	// Returns the energy value given deformed `x` and undeformed `X` state
	virtual double getEnergy(const VectorXd& x, const VectorXd& X) {
#ifdef DEBUG_OLD
      std::cout << "\n> DEBUG_OLD getEnergy - size of x (input) : " << x.rows() << " by " << x.cols() <<  std::endl;
      std::cout << "> DEBUG_OLD getEnergy - size of X (input) : " << X.rows() << " by " << X.cols() <<  std::endl;
#endif

		assert(x.rows() == X.rows());

		// Ex 1.2
		//Compute actual length
		Vector2d elong_Pos_0 = getNodePos(0, x);
		Vector2d elong_Pos_1 = getNodePos(1, x);
		double l = getLength(elong_Pos_0, elong_Pos_1);
	
		//Compute rest length
		Vector2d rest_Pos_0 = getNodePos(0, X);
		Vector2d rest_Pos_1 = getNodePos(1, X);
		double L = getLength(rest_Pos_0, rest_Pos_1);

		// Energy = 0.5 * k * ((l - L)/L)²
		return 0.5 * k *L* ((l-L)/L) * ((l-L)/L);

#ifdef DEBUG_OLD
      std::cout << "\n> DEBUG_OLD getEnergy - l (inside) : " << l <<  std::endl;
      std::cout << "> DEBUG_OLD getEnergy - L (inside) : " << L <<  std::endl;
      std::cout << "> DEBUG_OLD getEnergy - elongation (inside) : " << elongation <<  std::endl;
      std::cout << "> DEBUG_OLD getEnergy - energy (output) : " << 0.5*k*elongation*elongation <<  std::endl;
#endif
		//Compute the energy
		// 0.5 * k * (l - Lrepos)²

		// Some notes:
		// `x` and `X` contain the current and rest positions of all
		// nodes. You can extract the position of e.g. node 0 like this:
		// Vector2d x1 = getVertex(0, x);
		// or to get the rest position of node 0:
		// Vector X1 = getVertex(0, X);
		// The spring stiffness is stored in the variable `k`.
	}


/*
Wolfram : 
sqrt((x_i-x_j)²+(y_i-y_j)²)
derivate -k* ( (l(x,y,z,w)/L)-1) * u(x,y,z,w)
derivate -k* ( (l(x,y,z,w)/L)-1) * (x-y)/l(x,y,z,w)
derivate -k* ( (l(x,y,z,w)/L)-1) * (z-w)/l(x,y,z,w)

derivate 0.5*k* ( (l(x,y,z,w)/L)-1)²
Hessian matrix 0.5*k* ( (l(x,y,z,w)/L)-1)²

derivate -k*((l(w,x,y,z)/L)-1)*u(w,x,y,z)
derivate (w-z)/sqrt((x-y)²+(w-z)²)

derivate -k*((    sqrt((w-x)²+(y-z)²)/L)-1)  *   ([(w-x)]/sqrt((w-x)²+(y-z)²)) according to w,x,y,z
derivate -k*((    sqrt((w-x)²+(y-z)²)/L)-1)  *   ([(y-z)]/sqrt((w-x)²+(y-z)²)) according to w,x,y,z

derivate 0.5*k*( (sqrt((w-x)^2+(y-z)^2)/L -1)² according to w,x,y,z
Hessian matrix 0.5*k*( (sqrt((w-x)^2+(y-z)^2)/L -1)²
*/
	// Adds the gradient to `grad` given deformed `x` and undeformed `X` state
	virtual void addEnergyGradientTo(const VectorXd& x, const VectorXd& X, VectorXd& grad) {
#ifdef DEBUG_OLD
      std::cout << "\n> DEBUG_OLD addEnergyGradientTo - size of x (input) : " << x.rows() << " by " << x.cols() <<  std::endl;
      std::cout << "> DEBUG_OLD addEnergyGradientTo - size of X (input) : " << X.rows() << " by " << X.cols() <<  std::endl;
      std::cout << "> DEBUG_OLD addEnergyGradientTo - size of grad (input) : " << grad.rows() << " by " << grad.cols() <<  std::endl;
#endif

		assert(x.rows() == X.rows());

		// Ex 1.2
		// Task: Given `x` and `X`, add the gradient of the spring energy to `grad`.

		//We consider the gradient of E as : f = -dE
		// f = - k * ((l/L) - 1) * u
		//Compute rest length
		Vector2d rest_Pos_A = getNodePos(0, X); 
		Vector2d rest_Pos_B = getNodePos(1, X); 
		double L_norm = getLength(rest_Pos_A, rest_Pos_B);

		//Compute actual length
		Vector2d elong_Pos_A = getNodePos(0, x); 
		Vector2d elong_Pos_B = getNodePos(1, x);
		Vector2d l_vec = elong_Pos_A - elong_Pos_B;
		double l_norm = getLength(elong_Pos_A,elong_Pos_B);

		Vector2d u = l_vec.normalized();

		double cste = k * ((l_norm-L_norm)/L_norm);
		assert(std::isfinite(cste));

		grad[2*getNodeIndex(0)] += cste * u[0];
		grad[2*getNodeIndex(0)+1] += cste * u[1];

		grad[2*getNodeIndex(1)] += -cste * u[0];
		grad[2*getNodeIndex(1)+1] += -cste * u[1];

#ifdef DEBUG
		std::cout << "Bloc Gradient : " << std::endl;
		std::cout << "cste * u[0]  \t" << cste * u[0] << " \t -cste * u[0] " << -cste * u[0] <<  std::endl;
		std::cout << "cste * u[1]  \t" << cste * u[1] << " \t -cste * u[1] " << -cste * u[1] <<  std::endl;
		numericalDiffGrad(x,X);
		std::cout << " ----------------- " << std::endl;
#endif

		// Again, you can extract the position of e.g. node 0 like this:
		// Vector2d x1 = getVertex(0, x);
		// and the spring stiffness is stored in `k`.

		// Remember that `grad` is a vector of size 2*N, where N is the total
		// number of nodes in the system. Make sure you are writing to the
		// correct location in `grad`. To get the global index of node 0 of
		// this spring, use this function:
		// int globalIndex0 = getNodeIndex(0);
		// or for node 1
		// int globalIndex1 = getNodeIndex(1);
	}

	virtual void numericalDiffGrad(const VectorXd& x, const VectorXd& X){
		const double CSTE_FACTOR = 1;

		//Get the machine precision
		double eps = std::numeric_limits<double>::epsilon()*CSTE_FACTOR;
		//double eps = 1e-8;
		
		// We have a current Energy value (for actual position) and we want to know "how much" it variates thanks to x component and thanks to y component.
		// We move on a very low scale the position on one axis, and compute th new energy of the system. The variation is the gradient
		// We consider the elongation with NodeIndex1 Fixed and nodeIndex2 movable

		double actual_f_value = 0; //this->getEnergy(x,X);

		//Compute rest length
		Vector2d rest_Pos_A = getNodePos(0, X); 
		Vector2d rest_Pos_B = getNodePos(1, X); 
		Vector2d L_vec = rest_Pos_A - rest_Pos_B;
		double L_norm = getLength(rest_Pos_A, rest_Pos_B);

		//Compute actual length
		Vector2d elong_Pos_A = getNodePos(0, x); 
		Vector2d elong_Pos_B = getNodePos(1, x);
		Vector2d l_vec = elong_Pos_A - elong_Pos_B;
		double l_norm = getLength(elong_Pos_A,elong_Pos_B);

		VectorXd elong_Pos_A_xPlus(2);
		VectorXd elong_Pos_A_xMinus(2);
		VectorXd elong_Pos_A_yPlus(2);
		VectorXd elong_Pos_A_yMinus(2);

		VectorXd elong_Pos_B_xPlus(2);
		VectorXd elong_Pos_B_xMinus(2);
		VectorXd elong_Pos_B_yPlus(2);
		VectorXd elong_Pos_B_yMinus(2);

		elong_Pos_A_xPlus = elong_Pos_A;
		elong_Pos_A_xMinus = elong_Pos_A;
		elong_Pos_A_yPlus = elong_Pos_A;
		elong_Pos_A_yMinus = elong_Pos_A;

		elong_Pos_B_xPlus = elong_Pos_A;
		elong_Pos_B_xMinus = elong_Pos_A;
		elong_Pos_B_yPlus = elong_Pos_A;
		elong_Pos_B_yMinus = elong_Pos_A;

		//Populate the vector, to calculate Energy(x+eps,y), Energy(x-eps,Y) ...
		elong_Pos_A_xPlus[0] += eps;
		elong_Pos_A_xMinus[0] -= eps;
		elong_Pos_A_yPlus[1] += eps;
		elong_Pos_A_yMinus[1] -= eps;

		elong_Pos_B_xPlus[0] += eps;
		elong_Pos_B_xMinus[0] -= eps;
		elong_Pos_B_yPlus[1] += eps;
		elong_Pos_B_yMinus[1] -= eps;

		//New springs
		Vector2d l_vec_A_xPlus = elong_Pos_A_xPlus - elong_Pos_B;
		Vector2d l_vec_A_xMinus = elong_Pos_A_xMinus - elong_Pos_B;
		Vector2d l_vec_A_yPlus = elong_Pos_A_yPlus - elong_Pos_B;
		Vector2d l_vec_A_yMinus = elong_Pos_A_yMinus - elong_Pos_B;
		Vector2d l_vec_B_xPlus = elong_Pos_A - elong_Pos_B_xPlus;
		Vector2d l_vec_B_xMinus = elong_Pos_A - elong_Pos_B_xMinus;
		Vector2d l_vec_B_yPlus = elong_Pos_A - elong_Pos_B_yPlus;
		Vector2d l_vec_B_yMinus = elong_Pos_A - elong_Pos_B_yMinus;
		//Vector2d gradL = computeGrad(k,l,L_norm,u);

		//Calculate the directionnal derivatives
		double Edxi = (computeEnergy(l_vec_A_xPlus,L_norm)-computeEnergy(l_vec_A_xMinus,L_norm))/(2*eps);
		double Edyi = (computeEnergy(l_vec_A_yPlus,L_norm)-computeEnergy(l_vec_A_yMinus,L_norm))/(2*eps);
		double Edxj = (computeEnergy(l_vec_B_xPlus,L_norm)-computeEnergy(l_vec_B_xMinus,L_norm))/(2*eps);
		double Edyj = (computeEnergy(l_vec_B_yPlus,L_norm)-computeEnergy(l_vec_B_yMinus,L_norm))/(2*eps);

	#ifdef DEBUG
			std::cout << "Bloc finite difference (grad E): " << std::endl;
			std::cout << "Edxi \t" << Edxi << "\t Edxj \t" << Edxj << std::endl;
			std::cout << "Edyi \t" << Edyi << "\t Edyj \t" << Edyj << std::endl;
	#endif

	}

	virtual void numericalDiffHess(const VectorXd& x, const VectorXd& X){
		const double CSTE_FACTOR = 4;

		//Get the machine precision
		double eps = std::numeric_limits<double>::epsilon()*CSTE_FACTOR;

		// We have a current Energy value (for actual position) and we want to know "how much" it variates thanks to x component and thanks to y component.
		// We move on a very low scale the position on one axis, and compute th new energy of the system. The variation is the gradient
		// We consider the elongation with NodeIndex1 Fixed and nodeIndex2 movable

		double actual_f_value = 0; //this->getEnergy(x,X);

		//Compute rest length
		Vector2d rest_Pos_A = getNodePos(0, X); 
		Vector2d rest_Pos_B = getNodePos(1, X); 
		Vector2d L_vec = rest_Pos_A - rest_Pos_B;
		double L_norm = getLength(rest_Pos_A, rest_Pos_B);

		//Compute actual length
		Vector2d elong_Pos_A = getNodePos(0, x); 
		Vector2d elong_Pos_B = getNodePos(1, x);
		Vector2d l_vec = elong_Pos_A - elong_Pos_B;
		double l_norm = getLength(elong_Pos_A,elong_Pos_B);

		VectorXd elong_Pos_A_xPlus(2);
		VectorXd elong_Pos_A_xMinus(2);
		VectorXd elong_Pos_A_yPlus(2);
		VectorXd elong_Pos_A_yMinus(2);

		VectorXd elong_Pos_B_xPlus(2);
		VectorXd elong_Pos_B_xMinus(2);
		VectorXd elong_Pos_B_yPlus(2);
		VectorXd elong_Pos_B_yMinus(2);

		elong_Pos_A_xPlus = elong_Pos_A;
		elong_Pos_A_xMinus = elong_Pos_A;
		elong_Pos_A_yPlus = elong_Pos_A;
		elong_Pos_A_yMinus = elong_Pos_A;

		elong_Pos_B_xPlus = elong_Pos_A;
		elong_Pos_B_xMinus = elong_Pos_A;
		elong_Pos_B_yPlus = elong_Pos_A;
		elong_Pos_B_yMinus = elong_Pos_A;

		//Populate the vector, to calculate Energy(x+eps,y), Energy(x-eps,Y) ...
		elong_Pos_A_xPlus[0] += eps;
		elong_Pos_A_xMinus[0] -= eps;
		elong_Pos_A_yPlus[1] += eps;
		elong_Pos_A_yMinus[1] -= eps;

		elong_Pos_B_xPlus[0] += eps;
		elong_Pos_B_xMinus[0] -= eps;
		elong_Pos_B_yPlus[1] += eps;
		elong_Pos_B_yMinus[1] -= eps;

		//New springs
		Vector2d l_vec_A_xPlus = elong_Pos_A_xPlus - elong_Pos_B;
		Vector2d l_vec_A_xMinus = elong_Pos_A_xMinus - elong_Pos_B;
		Vector2d l_vec_A_yPlus = elong_Pos_A_yPlus - elong_Pos_B;
		Vector2d l_vec_A_yMinus = elong_Pos_A_yMinus - elong_Pos_B;
		Vector2d l_vec_B_xPlus = elong_Pos_A - elong_Pos_B_xPlus;
		Vector2d l_vec_B_xMinus = elong_Pos_A - elong_Pos_B_xMinus;
		Vector2d l_vec_B_yPlus = elong_Pos_A - elong_Pos_B_yPlus;
		Vector2d l_vec_B_yMinus = elong_Pos_A - elong_Pos_B_yMinus;

		Vector2d u = l_vec.normalized();
		Vector2d u_A_xPlus = l_vec_A_xPlus.normalized();
		Vector2d u_A_xMinus = l_vec_A_xMinus.normalized();
		Vector2d u_A_yPlus = l_vec_A_yPlus.normalized();
		Vector2d u_A_yMinus = l_vec_A_yMinus.normalized();
		Vector2d u_B_xPlus = l_vec_B_xPlus.normalized();
		Vector2d u_B_xMinus = l_vec_B_xMinus.normalized();
		Vector2d u_B_yPlus = l_vec_B_yPlus.normalized();
		Vector2d u_B_yMinus = l_vec_B_yMinus.normalized();

		//Vector2d gradL = computeGrad(k,l,L_norm,u);

		//Calculate the directionnal derivatives
		Vector2d Fdxi = (computeGrad(k,l_vec_A_xPlus,L_norm,u_A_xPlus)-computeGrad(k,l_vec_A_xMinus,L_norm,u_A_xMinus))/(2*eps);
		Vector2d Fdyi = (computeGrad(k,l_vec_A_yPlus,L_norm,u_A_yPlus)-computeGrad(k,l_vec_A_yMinus,L_norm,u_A_yMinus))/(2*eps);
		Vector2d Fdxj = (computeGrad(k,l_vec_B_xPlus,L_norm,u_B_xPlus)-computeGrad(k,l_vec_B_xMinus,L_norm,u_B_xMinus))/(2*eps);
		Vector2d Fdyj = (computeGrad(k,l_vec_B_yPlus,L_norm,u_B_yPlus)-computeGrad(k,l_vec_B_yMinus,L_norm,u_B_yMinus))/(2*eps);

	#ifdef DEBUG
			std::cout << "Bloc finite difference : " << std::endl;
			std::cout << "Fdxi[0] \t" << Fdxi[0] << "\t Fdyi[0] \t" << Fdyi[0] << std::endl;
			std::cout << "Fdxi[1] \t" << Fdxi[1] << "\t Fdyi[1] \t" << Fdyi[1] << std::endl;
			std::cout << "Fdxj[0] \t" << Fdxj[0] << "\t Fdyj[0] \t" << Fdyj[0] << std::endl;
			std::cout << "Fdxj[1] \t" << Fdxj[1] << "\t Fdyj[1] \t" << Fdyj[1] << std::endl;
	#endif
	}

	// Adds the hessian entries to `hesEntries` given deformed `x` and undeformed `X` state
	virtual void addEnergyHessianTo(const VectorXd& x, const VectorXd& X, std::vector<Tripletd>& hesEntries) {

		// Ex 1.4
		// Task: Given `x` and `X`, add the hessian of the spring energy to `hesEntries`.
		
		//Compute rest length
		Vector2d rest_Pos_A = getNodePos(0, X); 
		Vector2d rest_Pos_B = getNodePos(1, X); 
		Vector2d L_vec = rest_Pos_A - rest_Pos_B;
		double L_norm = getLength(rest_Pos_A, rest_Pos_B);

		//Compute actual length
		Vector2d elong_Pos_A = getNodePos(0, x); 
		Vector2d elong_Pos_B = getNodePos(1, x);
		Vector2d l_vec = elong_Pos_A - elong_Pos_B;
		double l_norm = getLength(elong_Pos_A,elong_Pos_B);

		Vector2d u = l_vec.normalized();

		double ga = k/(L_norm*L_norm*l_norm*l_norm); //k/(L²*l²)
		double be = -(k*((l_norm/L_norm)-1))/(L_norm*l_norm*l_norm*l_norm); //-k*((l/L)-1)/(L*l^3)
		double al = (k*((l_norm/L_norm)-1))/(L_norm*l_norm); //k* (l/L -1)/(L*l)

		double A = l_vec[0]*l_vec[0]*ga + l_vec[0]*l_vec[0]*be + al; 
		double B = l_vec[0]*l_vec[1]*ga + l_vec[0]*l_vec[1]*be; 
		//double C = ; 
		double D = l_vec[1]*l_vec[1]*ga + l_vec[1]*l_vec[1]*be + al; 

		//Fix due to the error in the energy
		A *= L_norm;
		B *= L_norm;
		D *= L_norm;

		//Create hessian
		int nodeAIndex = 2*getNodeIndex(0);
		int nodeBIndex = 2*getNodeIndex(1);
		hesEntries.push_back(Tripletd(nodeAIndex, nodeAIndex, 			A )); 		// Haut gauche 	// Ligne haut, Exi, dxi
		hesEntries.push_back(Tripletd(nodeAIndex, nodeAIndex + 1,		B )); 	// Haut droite 	// Ligne haut, Exi, dyi
		hesEntries.push_back(Tripletd(nodeAIndex+1, nodeAIndex, 		B ));		// Bas gauche 	// Ligne bas, Exi, dxi
		hesEntries.push_back(Tripletd(nodeAIndex+1, nodeAIndex +1, 		D ));		// Bas droite 	// Ligne bas, Exi, dyi

		hesEntries.push_back(Tripletd(nodeAIndex, nodeBIndex, 			- A )); 	// Haut gauche 	// Ligne haut, Exi, dxj
		hesEntries.push_back(Tripletd(nodeAIndex, nodeBIndex + 1, 		-B )); 		// Haut droite 	// Ligne haut, Exi, dyj
		hesEntries.push_back(Tripletd(nodeAIndex+1, nodeBIndex,  		-B ));		// Bas gauche 	// Ligne bas, Exi, dxj
		hesEntries.push_back(Tripletd(nodeAIndex+1, nodeBIndex +1, 		- D ));		// Bas droite 	// Ligne bas, Exi, dyj

		hesEntries.push_back(Tripletd(nodeBIndex, nodeAIndex, 			- A )); 	// Haut gauche 	// Ligne haut, Exj,  dxi
		hesEntries.push_back(Tripletd(nodeBIndex, nodeAIndex + 1, 		-B )); 		// Haut droite 	// Ligne haut, Exj, dyi
		hesEntries.push_back(Tripletd(nodeBIndex+1, nodeAIndex,  		-B ));		// Bas gauche 	// Ligne bas, Exj, dxi
		hesEntries.push_back(Tripletd(nodeBIndex+1, nodeAIndex +1, 		- D ));		// Bas droite 	// Ligne bas, Exj, dyi

		hesEntries.push_back(Tripletd(nodeBIndex, nodeBIndex, 			A )); 		// Haut gauche 	// Ligne haut, Exj, dxj
		hesEntries.push_back(Tripletd(nodeBIndex, nodeBIndex + 1, 		 B )); 	// Haut droite 	// Ligne haut, Exj, dyj
		hesEntries.push_back(Tripletd(nodeBIndex+1, nodeBIndex,  		 B ));		// Bas gauche 	// Ligne bas, Exj, dxj
		hesEntries.push_back(Tripletd(nodeBIndex+1, nodeBIndex +1, 		D ));		// Bas droite 	// Ligne bas, Exj, dyj

#ifdef DEBUG
		std::cout << "Bloc Hessian : " << std::endl;
		std::cout << A << " \t" << -A << " \t" << B << " \t" << -B <<  std::endl;
		std::cout << -A << " \t" << A << " \t" << -B << " \t" << B <<  std::endl;
		std::cout << B << " \t" << -B << " \t" << D << " \t" << -D <<  std::endl;
		std::cout << -B << " \t" << B << " \t" << -D << " \t" << D <<  std::endl;
		numericalDiffHess(x,X);
#endif

	}

protected:
	// the collection of nodes that define the triangle element
	std::array<int, 2> nodeIndices;
	// spring stiffness
	double k = 20.0;
};
