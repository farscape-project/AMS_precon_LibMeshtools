# AMS_precon_LibMeshtools
An attempt at implementing AMS preconditioning within LibMesh 

# Aims
Implement the AMS-preconditioner from Hypre:
 - Implement the G-Discrete gradient operator in LibMesh
 - Implement a test case for comparison (curl-curl problem)


# Problems with LibMesh
LibMesh is a finite element library that uses element classes
to solve finite element problems for the nodes. As such it does not
store mesh topological information for other element sub-entities e.g.
faces, edges and volumes. It instead stores them relative to mappings
of nodes. Which means to gain topological information about these
sub-entities on-the-fly calculation has to be carried out to gain access to
it.


Exotic elements such as Nedlec, Raviart thomas, mixed-elements
which explicitly solve for values at these sub-entities cannot store 
information easily. This problem can be resolved partially 
by using lagrange multipliers to represent the sub-entities and adding
additional equations, however when it comes to coarsening of the sub-entity spaces
the topological inforation may be relevant.

# Problem
I don't aim to majorly change the library, its too mature and I am too tired
i'm just going to generate an on-the-fly G-discrete-operator matrix and go from there.

