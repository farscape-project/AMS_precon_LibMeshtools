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
// Unique edges version
//
class G_operator
{
  private:
    // A map containing all the element edges
    // and mapping them to a pair of nodes (the
    // endpoints of the edge)
    std::map<int,std::pair<unsigned int, unsigned int>> edge_map;

    // Total number of local edges
    unsigned int ntot_edges_local = 0;

    // Total number of global edges
    unsigned int ntot_edges_global = 0;

    //Hypre matrix objects
    HYPRE_IJMatrix     par_G_ij;
    HYPRE_ParCSRMatrix par_G;
    unsigned int ilower, iupper; //Edge lower and upper bounds
    unsigned int jlower, jupper; //Node lower and upper bounds

    //MPI processor communication stuff
    //MPI_Comm  comm;
    unsigned int nmessages=0;
    bool is_parallel=false;
    int ier, nprocs, procID;
    std::vector<int> ProcNeighbors;
    std::vector<unsigned int> ProcEdgeSize;

    //Function that takes the local edge number
    //and makes it global
    unsigned int local_to_global_edge(unsigned int local_edge);

    // Forms a list of the local process neighbors using the
	// nodal partitioning tables
    void Map_Proc_Neighbors(EquationSystems & es);

    // Makes the edge map and calculates number
    // of local and global edges
    void Make_Edge_Map(EquationSystems & es);

    // Called within make edgemap if the system is parallel
    // removes all remote edges and assignns unique edges
	// a unique proc ID
    void prune_Remote_Duplicate_Edges();

    // Sizes up the G-operator matrix for the PETSc-hypre 
    // interface
    void Size_G_Operator();

    // Sets the G-operator matrix using the PETSc-hypre 
    // interface
    void Set_G_Operator();

    //Cantor counting function
    int Cantors_Counter(unsigned int I , unsigned int J);

    //Inverse Cantor counting function
    std::pair<unsigned int, unsigned int> Cantors_CounterInv(int K);

    //Finds the difference vector between a pair of vectors
    template<typename T>
    void VectorDifference(std::vector<T> c, std::vector<T> a, std::vector<T> b){
      c.clear()
      if(a.size() == b.size()) for(unsigned int I=0; I<a.size(); I++) c.push_back(a[I] - b[I]);
      if(a.size() != b.size()) std::cout << "Error vectors don't agree in size"
    };

    //sign of inner product of a pair of vectors
    template<typename T>
    T InnerProductSign(std::vector<T> a, std::vector<T> b){
      if(a.size() != b.size()){ //making sure the vector sizes agree
        std::cout << "Error vectors don't agree in size"
        return T(0);
	  }
      T ab = T(0), signab = T(0);
      for(unsigned int I=0; I<a.size(); I++) ab = ab + a[I]*b[I];
      if(ab >  T(0) ) signab = T( 1);
      if(ab <= T(0) ) signab = T(-1);
      return signab;
    };

  public:
    G_operator(EquationSystems & es);
}

