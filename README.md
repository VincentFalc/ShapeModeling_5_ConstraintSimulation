---
author: VF
layout: post-full
type: image
featimg: GradientDescent.gif
title: Shape Modeling Forces simulation
tags: [Shape Modeling, IGL, Eigen, Simulation, Forces]
category: [Shape Modeling]
---

##### TL;DR
This project was linked to a lecture about Shape Modeling. This project presents a method to obtain the stable position of a mesh constrained under forces. The simulation is mainly based on Springs.
<br/>
The main idea is to compare the mesh with a Spring's network. Then, the problem is to minimize the global stress (which's the sum of the individual Spring stress) of the mesh, by an optimization method. Here again we are talking about reducing the energy of the mesh.
<br/>
We have choices on the Spring modelisation and the optimization method. We will go through first order optimization (Gradient descent) and second order optimization (Newton's method).
<br/>

##### Table of Contents  
[1. Mass-spring simulation](#Mass)  <br/>
[1.1. Gradient Descent and Line Search](#Gradient)  <br/>
[1.2. Spring simulation with gradient descent](#simulation)  <br/>
[1.3. Newton’s method](#Newton)  <br/>
[1.4. Spring simulation with Newton’s method](#method)  <br/>
[2. FEM Simulation](#FEM)  <br/>
[2.1. FEM simulation with gradient descent](#Graddescent)  <br/>
[2.2. FEM simulation with Newton’s method](#Newdescent)  <br/>

<a name="Mass"/>

#### 1. Mass-spring simulation

<a name="Gradient"/>

#### 1.1. Gradient Descent and Line Search

We can validate our gradient descent method on small and easy to debug problem.
So for this first part, we're going to validate the method on the Rosenbrock function
[-> Function](https://en.wikipedia.org/wiki/Rosenbrock_function)
<br/>
We will write several function : 
* computeValue(const VectorXd& x) will return f(x), with f the Rosenbrock function
* addGradientTo(VectorXd& grad, const VectorXd& x) will add to the gradient vector (first paramter) the computed gradient of the Rosenbrock function (analytical)
* doLineSearch will performs a line search of a better (lower value) solution, in a particular direction (given by dx) of research
<br/>
###### Explanation
We do compute the value as a simple function evaluation.

We do calculate the gradient of the function anaticaly and then compute it at one specific point. We get the two directional derivative of this function at that point.

<p align="center">
  <img height="200" src="./results/Gradient.png">
</p>

We do compute the search direction as the opposite of the gradient.
We do a line search as the current position, added of a 'dx' in the search direction, that we halves if the evaluation of the function at this position is higher than the current point - meaning that it makes no sense to "jump" there.
So, regarding the first value of Dx, we will consider 1*Dx,0.5*Dx,0.25*Dx, 0.125*Dx ... until we get a "better" value of the total Energy.

###### Usage
```C++
    virtual double computeValue(const VectorXd& x) {
		// Ex 1.1
		double xV = x[0];
		double yV = x[1];
		return (a-xV) * (a-xV) + b * (yV-(xV*xV))*(yV-(xV*xV));
	}

    (...)

    virtual void addGradientTo(VectorXd& grad, const VectorXd& x) {

		// Ex 1.1
		double xV = x[0];
		double yV = x[1];

		//Calculate the directionnal derivatives
		double dfdx = 2 * (-a + xV + 2 * b * xV * (xV*xV - yV));
		double dfdy = 2 * b * (yV - xV*xV);

		//Store it
		grad[0] += dfdx; //+ grad[0];
		grad[1] += dfdy; //+ grad[1]; /// ? Adds ? Or replace ? Or appends ?
    }

    (...)

    virtual void computeSearchDirection(ObjectiveFunction *function, const VectorXd &x, VectorXd &dx)
	{
		// Ex. 1.1
		VectorXd grad;
		grad.resize(x.rows(), 1);
		function->addGradientTo(grad, x);

		//Store it
		dx[0] = -grad[0];
		dx[1] = -grad[1];
		// Note : no normalization, otherwise we loose information ?
	}

    (...)

	virtual void doLineSearch(ObjectiveFunction *function, const VectorXd &dx, VectorXd &xi)
	{

		// Ex. 1.1
        (...)

		for (int i = 0; i < maxLineSearchIterations && !isInferior; i++)
		{

			if (function->computeValue(tmpx) > actualValue)
			{ //We have " a higher point "
				tmpx = xi + dx*currentAlpha;
				currentAlpha *= beta;
			}
			else
			{ // We have a lower point, we stop here
				isInferior = true; //Not really necessary ...
				break;
			}
		}

		//Storage
		xi = tmpx;
	}
    (...)

```

###### Result

The minimum of the Rosenbrock functin is (1,1) and we see the gradient descent find this value, but with a high number of iterations.

```
Gradient descent:
converged:  yes
iterations: 9750
min. x:     1.00001 1.00002
min. value: 6.17752e-11
```

__________________________________

<a name="simulation"/>

#### 1.2. Spring simulation with gradient descent

Now that your basic function are functionnal, we can get to a harder problem : the spring simulation of the mesh.
<br/>

Similarly, we will write several function : 
* getEnergy(const VectorXd& x, const VectorXd& X) will spring energy given the current list of state x and list of rest state X of every spring.
* addEnergyGradientTo(const VectorXd& x, const VectorXd& X, VectorXd&
grad), which add to the gradient vector (last paramter) the computed gradient of the spring energy given states x and rest states X.
* doLineSearch will be the same as previously.
<br/>

###### Explanation
We calculate the energy of a Spring as the expression provided in the PDF.

We calculate the gradient of the energy as : f = -dE. So we calculate the opposite of the Force function, as provided in the PDF, at each point (2 ends) of the current Spring.
The calculated value is stored in the gradient Vector, at the place of each node (x and y components).

###### Usage
```C++
	virtual double getEnergy(const VectorXd& x, const VectorXd& X) {
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
		return 0.5 * k * L * (l-L)/L * (l-L)/L;
	}

	virtual void addEnergyGradientTo(const VectorXd& x, const VectorXd& X, VectorXd& grad) {
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

	}

```

###### Result
<div style="text-align:center">
<p align="center">
<img height="300" src="./results/GradientDescent.png">
<img height="300" src="./results/GradientDescent.gif">
<br/>
  Result of "Test" button - Gradient Descent - Spring - Default max stress (0.02)
</p>
</div>

The gradient descent give a minimal value, with a high number of iteration, that we will be able to compare with output of other methods.

```
Gradient Descent
total energy = 117.884
# iterations = 3468
min def at   = 7.80134e-08
max def at   = 0.00641764
```

__________________________________
<a name="Newton"/>

#### 1.3. Newton’s method

The first order can get to a solution, but need many iterations (especially in some particular situation, if the function is of degree 2 or more).
<br/>
Therefore, we can compute a second order solution. As for the first order, we begin by the Rosenbrock Function.
<br/>
We will write several function : 
* addHessianEntriesTo, which adds its second order derivatives as Eigen::Triplets.
* computeSearchDirection, which computes the search direction dx given an objective function. We use Eigen::SimplicialLDLT to solve the search Direction.

###### Explanation
We calculate analytically the hessian, and we compute it directly.

<p align="center">
  <img height="500" src="./results/hessian.jpg">
</p>

We then solve the system to get the dx value.

###### Usage
```C++
virtual void addHessianEntriesTo(std::vector<Tripletd> &hessianEntries, const VectorXd &x)
	{
		// Ex 1.2
		// write d^2f/dx^2 in `hessianEntries`
		double xV = x[0];
		double yV = x[1];

		//Fill the Hessian
		hessianEntries.push_back(Tripletd(0, 0, -4 * b * (yV - xV * xV) + 8 * b * xV * xV + 2));
		hessianEntries.push_back(Tripletd(0, 1, -4 * b * xV));
		hessianEntries.push_back(Tripletd(1, 0, -4 * b * xV));
		hessianEntries.push_back(Tripletd(1, 1, 2 * b));
	}

    (...)

	virtual void computeSearchDirection(ObjectiveFunction *function, const VectorXd &x, VectorXd &dx)
	{
    (...)
		// Ex 1.3
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    (...)

		function->addHessianEntriesTo(hessianEntries, x);
		hessianSparse.setFromTriplets(hessianEntries.begin(), hessianEntries.end());

		function->addGradientTo(grad, x);
		grad = -1 * grad;
		rightSide = grad.sparseView();

		solver.compute(hessianSparse);
		solved = solver.solve(rightSide);

		dx = MatrixXd(solved);
	}

    (...)

```

###### Result

We see, compared to the gradient descent, that Newton method give a correct result, but only in one iteration. 
Note : the 0.000001 of error for the gradient descent is due to the way the calculation is stopped. If we don't see any improvement (0.00001 or more difference) between two iterations, we stop the process. As the process converge very slowly, there is not more improvement (numerically) and so the process is stopped.
Whereas for the Newton's Method, in one jump, it goes right at the good position.

```
Newton's method:
converged:  yes
iterations: 4
min. x:     1 1
min. value: 2.91398e-13
```

__________________________________
<a name="method"/>

#### 1.4. Spring simulation with Newton’s method


Now that your basic function are functionnal for the second order method, we can get similarly improve it to the Spring simulation.
<br/>

Similarly, we will write several function : 
* addHessianEntriesTo(const VectorXd& x, const VectorXd& X, VectorXd&
grad), adds the second order derivatives of the spring energy to the Hessian (hesEntries)
<br/>

###### Explanation
We compute analytically the Hessian and implement it the same way. The disposition in the "global" hessian is dependent on which point we consider and which derivative of it we compute. 

The manual computation is verified automatically : 
<p align="center">
  <img height="700" src="./results/Hessian.png">
</p>

###### Usage
```C++
virtual void addHessianEntriesTo(std::vector<Tripletd> &hessianEntries, const VectorXd &x)
	{
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

    (...)

		//Create hessian
		int nodeAIndex = 2*getNodeIndex(0);
		int nodeBIndex = 2*getNodeIndex(1);
		hesEntries.push_back(Tripletd(nodeAIndex, nodeAIndex, 			A )); 	// Haut gauche 	// Ligne haut, Exi, dxi
		hesEntries.push_back(Tripletd(nodeAIndex, nodeAIndex + 1,		B )); 	// Haut droite 	// Ligne haut, Exi, dyi
		hesEntries.push_back(Tripletd(nodeAIndex+1, nodeAIndex, 		B ));	// Bas gauche 	// Ligne bas, Exi, dxi
		hesEntries.push_back(Tripletd(nodeAIndex+1, nodeAIndex +1, 		D ));	// Bas droite 	// Ligne bas, Exi, dyi

		hesEntries.push_back(Tripletd(nodeAIndex, nodeBIndex, 			- A )); // Haut gauche 	// Ligne haut, Exi, dxj
		hesEntries.push_back(Tripletd(nodeAIndex, nodeBIndex + 1, 		-B )); 	// Haut droite 	// Ligne haut, Exi, dyj
		hesEntries.push_back(Tripletd(nodeAIndex+1, nodeBIndex,  		-B ));	// Bas gauche 	// Ligne bas, Exi, dxj
		hesEntries.push_back(Tripletd(nodeAIndex+1, nodeBIndex +1, 		- D ));	// Bas droite 	// Ligne bas, Exi, dyj

		hesEntries.push_back(Tripletd(nodeBIndex, nodeAIndex, 			- A )); // Haut gauche 	// Ligne haut, Exj,  dxi
		hesEntries.push_back(Tripletd(nodeBIndex, nodeAIndex + 1, 		-B )); 	// Haut droite 	// Ligne haut, Exj, dyi
		hesEntries.push_back(Tripletd(nodeBIndex+1, nodeAIndex,  		-B ));	// Bas gauche 	// Ligne bas, Exj, dxi
		hesEntries.push_back(Tripletd(nodeBIndex+1, nodeAIndex +1, 		- D ));	// Bas droite 	// Ligne bas, Exj, dyi

		hesEntries.push_back(Tripletd(nodeBIndex, nodeBIndex, 			A )); 	// Haut gauche 	// Ligne haut, Exj, dxj
		hesEntries.push_back(Tripletd(nodeBIndex, nodeBIndex + 1, 		 B )); 	// Haut droite 	// Ligne haut, Exj, dyj
		hesEntries.push_back(Tripletd(nodeBIndex+1, nodeBIndex,  		 B ));	// Bas gauche 	// Ligne bas, Exj, dxj
		hesEntries.push_back(Tripletd(nodeBIndex+1, nodeBIndex +1, 		D ));	// Bas droite 	// Ligne bas, Exj, dyj

	}

```

###### Result
<div style="text-align:center">
<p align="center">
  <img height="300" src="./results/NewtonDescent.png">
  <img height="300" src="./results/NewtonDescent.gif">
</p>
</div>

Here again, we can see that's the second order converge much faster, to the same result.

```
total energy = 117.539
# iterations = 3
min def at   = 2.53686e-11
max def at   = 0.0172078
```

__________________________________
<a name="FEM"/>

### 2. FEM Simulation

<a name="Graddescent"/>

#### 2.1. FEM simulation with gradient descent

If the Spring modelization is easy enough, as previously, we can compute derivative analyticaly. In case of more complex Spring modelization - as Neo-Hookean Springs - we can't easily analyticaly calculate the derivative, as it is already of degree 2.
<br/>
In this case, we might use numerical derivation, or automatic differentiation.
We're going to see an implementation of the numerical derivation solution, with Finite Element Method.
<br/>
Some functions will be implemented : 
* getEnergy, which returns the Neo-Hookean deformation energy of one Spring.
* addEnergyGradientTo, which write the gradient of this deformation in the given array.

###### Explanation
We compute the energy from the energy density, by simply multiplicate the energy density by the area where this energy is applied. The expression of the Energy density is the Neo-Hookean one.
<br/>
We then compute the Gradient, with finite differences, as the difference of energy, normalized, if points involved in the structure are sligtly moved.

###### Usage
```C++

	virtual double getEnergy(const VectorXd &x, const VectorXd &X)
	{
		// Ex 2.1
		Vector2d xTMP[3];
		xTMP[0] = getNodePos(0, x);
		xTMP[1] = getNodePos(1, x);
		xTMP[2] = getNodePos(2, x);

		computeDeformationGradient(xTMP, dxdX);
		MatrixXd C = (dxdX.transpose() * dxdX);
		double logdetF = log(dxdX.determinant());

		double Energy = (shearModulus * 0.5 * (C.trace() - 2) - shearModulus * logdetF + bulkModulus * 0.5 * logdetF * logdetF);
		return restShapeArea * Energy;
	}

	virtual void addEnergyGradientTo(const VectorXd &x, const VectorXd &X, VectorXd &grad)
	{
		// Ex 2.1
		for (int j = 0; j < 2; j++)
		{
			for (int i = 0; i < 3; i++)
			{
				modified_x_A = x;
				modified_x_B = x;

				modified_x_A[2 * indexPoint[i] + j] += h;
				modified_x_B[2 * indexPoint[i] + j] -= h;
				double Etmp = (getEnergy(modified_x_A, X) - getEnergy(modified_x_B, X)) / (2 * h);
				grad[2 * indexPoint[i] + j] += Etmp;
			}
		}
```

###### Result
<div style="text-align:center">
<p align="center">
  <img height="300" src="./results/FEM_GradientDescent.png">
</p>
</div>

We get a solution, with a high number of iterations.

```
total energy = 115.291
# iterations = 6609
min def at   = 4.72342e-08
max def at   = 0.0114811
```

__________________________________
<a name="Newdescent"/>

#### 2.2. FEM simulation with Newton’s method

Here again, we can extend this derivation to the second order. 
We will need the function addEnergyHessianTo which will add the values of second order derivative to the provided matrix.

###### Explanation
We construct the Hessian with finite difference, as the variation of gradients, if we slighlty move 2 points of the initial structure.
As there is a lot (6*6) of coefficients to compute, we use iterative structures.
<br/>
The positions of each coefficient depends on the 2 points involved (that gives the X position and Y position of the 4*4 "block") and if we are considering X or Y coordinates (which give the position into that block, dxdx in upper left, dydy in bottom right).

###### Usage
```C++
virtual void addEnergyHessianTo(const VectorXd &x, const VectorXd &X, std::vector<Tripletd> &hesEntries)
{
	// Ex 2.2

	//Liste of necessary energies
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

```

###### Result
<div style="text-align:center">
<p align="center">
  <img height="300" src="./results/FEM_Newton.png">
</p>
</div>

Here again, we see a similar solution found (even slighlty better) with a much lower number of iterations.

```
total energy = 114.438
# iterations = 6
min def at   = 6.79233e-08
max def at   = 0.0321732
```


For more information, see : 
[-> Repository](https://github.com/VincentFalc/ShapeModeling_5_ConstraintSimulation)

Note that the code is not optimized by thoughtful laziness: this code is not made to be reused as is in other applications, and if this is the case, the code will be reviewed and optimized. Thanks for your understanding.
