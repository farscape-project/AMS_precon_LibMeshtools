

#include "G_operator_gen.hpp"




void EdgeMap::Make_Edge_Map(EquationSystems & es, const std::string & system_name)
{
  // Ignore unused parameter warnings when !LIBMESH_ENABLE_AMR.
  libmesh_ignore(es, system_name);

  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "System");

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();


  //=====
  // Form the initial edge-node-pair-map
  //=====
  ntot_edges_local = 0; //Just making sure its zero;
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
      edge_map[ntot_edges_local]= edge_key;
      ntot_edges_local++;
    }
  }

  //=====
  // Find the global number of edges 
  //=====
  ntot_edges_global = 0; //Just making sure its zero;
  MPI_Allreduce(&ntot_edges_local, &ntot_edges_global, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);


  //=====
  // Find the minimum and maximum global
  // node numbers of local process
  //=====
  for(it = edge_map.begin(); it != edge_map.end(); it++){
    
  }
}




void EdgeMap::Find_G_Operator(){
  for (it = edge_map.begin(); it != edge_map.end(); it++)
  {
    unsigned int k = it->first;               // Local Edge number
    unsigned int m = (it->second).first;      // first  Edge node (global numbering)
    unsigned int n = (it->second).second;     // second Edge node (global numbering)
    unsigned int l = local_to_global_edge(k); // Global Edge number



  }
};


// Sets the G-operator matrix using the PETSc-hypre 
// interface using the IJ matrix interface
void EdgeMap::Set_G_Operator(){

  for(unsigned int I=0; I<ntot_edges_local; I++){
    unsigned int K = edge_map[I].first;
    unsigned int L = edge_map[I].second;
  }



  HYPRE_IJMatrixCreate(comm, ilower, iupper, jlower, jupper, &par_G_ij);
};