# Sequential-Quadratic-Programming-method-Implementation-in-Matlab
A SQP algorithm implementation for solving nonlinear constrained optimization problems.

Summary of Steps for SQP Algorithm

1. Make a QP approximation to the original problem. For the first iteration, use a
Lagrangian Hessian equal to the identity matrix.

2. Solve for the optimum to the QP problem. As part of this solution, values for the
Lagrange multipliers are obtained.

3. Execute a simple line search by first stepping to the optimum of the QP problem. So the
initial step is ∆x, and xnew = xold + x. See if at this point a penalty function, composed of
the values of the objective and violated constraints, is reduced. If not, cut back the step
size until the penalty function is reduced. The penalty function is given by P = f + sum(λ*g), 
where the summation is done over the set of violated constraints, and
the absolute values of the constraints are taken. The Lagrange multipliers act as scaling
or weighting factors between the objective and violated constraints.

4. Evaluate the Lagrangian gradient at the new point. Calculate the difference in x and in
the Lagrangrian gradient, γ. Update the Lagrangian Hessian using the BFGS update.

5. Return to Step 1 until ∆x is sufficiently small. When ∆x approaches zero, the K-T
conditions for the original problem are satisfied.
