// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/sum_shell_matrix.h"
#include "libmesh/tensor_shell_matrix.h"
#include "libmesh/sparse_shell_matrix.h"

#include "libmesh/getpot.h"


// The definition of a geometric element
#include "libmesh/elem.h"
#ifndef HYPRE_AMS_INTERFACE_H
#define HYPRE_AMS_INTERFACE_H 

#include "libmesh/enum_solver_package.h"

#include <iostream>
#include <algorithm>
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>
#include <vector>
#include <pair>
#include <map>
#include <mpi.h>

//
// This class forms the AMS prconditioner components
// including the G-operator, direction vectors
// and pases it onto the Hypre/PETSc-hypre interface
//
class Hypre_AMS_Interface
{
  private:
    SupplementaryEntityIDs * _SupEiDs;

    // A map containing all the local edges
    // and mapping them to a pair of nodes (the
    // endpoints of the edge)
    std::map<int,std::pair<unsigned int, unsigned int>> edge_map;

    // Total number of local edges
    unsigned int ntot_edges_local = 0;

    // Total number of global edges
    unsigned int ntot_edges_global = 0;

    //Hypre linear algebra objects
	//Starting with the IJ object moving
    //into the Parallel CSR stored objects
    HYPRE_IJMatrix   par_G_ij;
    HYPRE_IJVector   x_ij_vec, y_ij_vec, z_ij_vec; //Coordinate vectors at vertices (IJ_vec)

    HYPRE_ParCSRMatrix par_G;                              //G-Operator CSR matrix
    HYPRE_ParVector    par_xcoord, par_ycoord, par_zcoord; //coordinates at vertices (CSR-vec)
    HYPRE_ParVector    par_xvec, par_yvec, par_zvec;       //Edge unit vectors (CSR-vec)

    unsigned int ilower, iupper; //Edge lower and upper bounds
    unsigned int jlower, jupper; //Node lower and upper bounds

    // Makes the edge map
    void Make_Edge_Map(EquationSystems & es, SupplementaryEntityIDs & SupEiDs);

    // Sets the G-operator matrix using the Hypre/PETSc-hypre
    // interface
    void Set_G_Operator();

  public:
    Hypre_AMS_Interface(EquationSystems & es, SupplementaryEntityIDs & SupEiDs);
}
#endif