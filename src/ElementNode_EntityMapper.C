class SupplementaryEntityIDs
{
    //These maps store a global contiguous numbering of element subentities
    //And map these to the global Node/Dof numbering system within LibMesh
    std::map<unsigned int, unsigned int> Vert_to_Global;
    std::map<unsigned int, unsigned int> Edge_to_Global;
    std::map<unsigned int, unsigned int> Face_to_Global;
    std::map<unsigned int, unsigned int> Volm_to_Global;
    unsigned int LocalEntitySizes[4] = {0,0,0,0};
    unsigned int LocalEntityStarts[4] = {0,0,0,0};

  public:

    template<typename S, typename T>
    void AddToMapIteratorIfUnique(std::map<S,T> EntityMap, S I, T J){
      if( EntityMap.find(J) == EntityMap.end() ){
        EntityMap[I] = J;
	    J++;
      }
    }

    //Form the entity Maps
    void FormEntityMaps(EquationSystems & es){
      // Get a constant reference to the mesh object.
      const MeshBase & mesh = es.get_mesh();

      //Find maps entities in global node numberings from local process
      //to a unique contiguous ID 
      unsigned int LProcID = mesh.processor_id();

      for (const auto & elem : mesh.active_local_element_ptr_range()){
        for(unsigned int I=0; I<elem->n_nodes(); I++){
          unsigned int nodeID = elem->node_id();

          Node & node = mesh.node_ref(nodeID);
          if( node.processor_id() != LProcID){
            if( elem->is_vertex(nodeID)   ) AddToMapIteratorIfUnique<unsigned int, unsigned int>(Vert_to_Global, nodeID, LocalEntitySizes[0]);
            if( elem->is_edge(nodeID)     ) AddToMapIteratorIfUnique<unsigned int, unsigned int>(Edge_to_Global, nodeID, LocalEntitySizes[1]);
            if( elem->is_face(nodeID)     ) AddToMapIteratorIfUnique<unsigned int, unsigned int>(Face_to_Global, nodeID, LocalEntitySizes[2]);
            if( elem->is_internal(nodeID) ) AddToMapIteratorIfUnique<unsigned int, unsigned int>(Volm_to_Global, nodeID, LocalEntitySizes[3]);
          }
        }
      }
    };
	//Nothing Interesting
};
