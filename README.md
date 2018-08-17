# Solute drag in Grain growth
 
This is a one dimensional C code, which can, at present simulate steady state one dimensional grain growth. 

CHANGELOG

06/07/18

The current code can simulate solute drag for an additional driving force in 1D and 2D considering solute drag. For more information on the model implemented, see Kim and Park (2008).  Memory error has been resolved and profile is almost equal to the Cahn's model (qualitatively). It now remains to study the variation of solute drag due to movement of the interface.

16/08/18

* Model correctly handles grain growth problem along with solute segregation in a single phase material.
* Constant diffusivity
* Slow in process
* Have to modify the gradient terms for multiphase field model
* Boundary conditions says drichlet but is actually neumann boundary condition
* Isothermal temperature

