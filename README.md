# CFD

## Introduction

This code is an implementation of the Finite Volume Method to numerically solve the Steady State 2D COnvection Diffusion equation. The code is written in MATLAB because of it's simplicity in protoyping numerical methods. However, it will be translated into Pyhton eventually. Currently the code only supports square domains and cartesian coordinates. There are other limitations that are outlined below.

## The Physics

Engineers make simulations like the one shown above by numerically solving the 2D Diffusion Equation over any given domain. For those who aren't familiar with the latter, the 2D Diffusion Equation is a Parabolic Partial Differential Equation (PDE), and in most cases it's numerically solved using either the Finite Volume Method (FVM) or the Finite Element Method (FEM). For a steady-state case (results do not vary with time) and no heat generation (there is no internal heat source in the domain) the 2D Diffusion Equation is defined as shown below:

$$ \frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} = 0 $$

## The Numerical Approach: Finite Volume Method

There are analytical solutions for the 2D Diffusion equation in simple cases such as the square domain of interest in this project, but analytically solving this PDE quickly becomes unfeasible with just slight increases in complexity, hence the need for numerical methods. Both the FEM and FVM approaches are based on the same philosophy: take the original continuous domain and break it up into small cells where the complicated PDE can be approximated by an algebraic equation in each one.

$$ T_{i+1,j} + T_{i-1,j} + T_{i,j+1} + T_{i,j-1} - 4T_{i,j} = 0 $$

The i and j subscripts refer to the temperature in a particular X and Y location of the mesh, i'e.: a particular cell.
