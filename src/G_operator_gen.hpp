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





void Make_Unique_Global_Edges(EquationSystems & es,
               const std::string & system_name)
{
  // Ignore unused parameter warnings when !LIBMESH_ENABLE_AMR.
  libmesh_ignore(es, system_name);

  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "System");

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();



  std::map<std::pair<unsigned int, unsigned int>, int> edge_keyMap
  edge_keyMap edge_map;
  unsigned int ntot_edges = 0;  


  //=====
  // Form the initial edge-node-pair-map
  //=====
  for (const auto & elem : mesh.active_local_element_ptr_range())
  {
    const unsigned int nedges = elem->n_edges();

    for(unsigned int I=0; I<nedges; I++)
    {
      //Find global nodeIDs of the edge endpoints
      unsigned int m = elem->node_id(elem->edge_nodes_map[I][0]);
      unsigned int n = elem->node_id(elem->edge_nodes_map[I][1]);

      //Form a unique pair that defines the edge
      // (inline with Nuno's nedlec numbering system)
      std::pair<unsigned int, unsigned int> edge_key;
      if(m < n) edge_key = make_pair(m,n);
      if(n < m) edge_key = make_pair(n,m);

      //Check whether the edge key exists within the
	  //edge-map
      if(edge_map.find(edge_key) == edge_map.end()){
        //If it does not exist add it to the end
		//otherwise do nothing
        edge_map[edge_key]=ntot_edges;
        ntot_edges++;
	  }
    }
  }


  //=====
  // Check for processor-shared nodes
  //=====


  //=====
  // Find edges containing those nodes
  //=====




  //=====
  // Remove duplicate entries
  //=====

}

void Find_G_Operator(){





};