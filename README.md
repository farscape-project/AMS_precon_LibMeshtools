# AMS_precon_LibMeshtools
An attempt at implementing AMS preconditioning within LibMesh 

# Aims
..Implement the AMS Hypre 



# Problems with LibMesh
LibMesh is a finite element library that uses element classes
to solve finite element problems for the nodes. As such it does not
store mesh topological information for other element sub-entities e.g.
faces, edges and volumes. It instead stores them relative to mappings
of nodes. Which means to gain topological information about these
sub-entities on-the-fly calculation has to be carried out to gain access to
it.
using 
and more exotic elements such as Nedlec, Raviart thomas, mixed-elements
which explicitly solve for values at these sub-entities cannot store 
information easily. And this is even harder
Some examples use lagrange multipliers and 
