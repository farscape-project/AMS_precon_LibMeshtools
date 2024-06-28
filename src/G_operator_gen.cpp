

#include "G_operator_gen.hpp"




void EdgeMap::Make_Unique_local_Edges(EquationSystems & es, const std::string & system_name)
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
        edge_map[edge_key]= ntot_edges; //Naieve numbering (needs updating)
        ntot_edges++;
	  }
    }
  }


  //=====
  // Find the global number of edges 
  //=====


  //=====
  // Update the edge numbers
  //=====
}




void EdgeMap::Find_G_Operator(){
  int k = 0;
  for (it = edge_map.begin(); it != edge_map.end(); it++)
  {
    k++;
    std::cout << (it->first).first;
              << it->second   // string's value 
              << std::endl;
}


};