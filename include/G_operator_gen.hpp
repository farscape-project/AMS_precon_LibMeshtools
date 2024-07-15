// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/sum_shell_matrix.h"
#include "libmesh/tensor_shell_matrix.h"
#include "libmesh/sparse_shell_matrix.h"
#include "libmesh/mesh_refinement.h"

#include "libmesh/getpot.h"


// The definition of a geometric element
#include "libmesh/elem.h"
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
// Non-unique edges version
//
class EdgeMap
{
  private:
    // A map containing all the element edges
    // and mapping them to a pair of nodes (the
    // endpoints of the edge)
    std::map<int,std::pair<unsigned int, unsigned int>> edge_map;

    // Iterator for the edge-map
    std::map<int,std::pair<unsigned int, unsigned int>>::iterator it;

    // Total number of local edges
    unsigned int ntot_edges_local = 0;

    // Total number of global edges
    unsigned int ntot_edges_global = 0;  

    //Hypre matrix objects
    HYPRE_IJMatrix     par_G_ij;
    HYPRE_ParCSRMatrix par_G;


    //Function that takes the local edge number
    //and makes it global
    unsigned int local_to_global_edge(unsigned int local_edge);

  public:
    // Makes the edge map and calculates number
    // of local and global edges (these are non-unique)
    void Make_Edge_Map(EquationSystems & es, const std::string & system_name);


    // Finds the necessary information needed to caculate
    // the G-operator matrix
    void Find_G_Operator();

    // Sets the G-operator matrix using the PETSc-hypre 
    // interface
    void Set_G_Operator();
}

