class SupplementaryEntityIDs
{
    //These maps store a global contiguous numbering of element subentities
    //And map these to the global Node/Dof numbering system within LibMesh
    std::map<unsigned int, unsigned int> Global_to_LVert;
    std::map<unsigned int, unsigned int> Global_to_LEdge;
    std::map<unsigned int, unsigned int> Global_to_LFace;
    std::map<unsigned int, unsigned int> Global_to_LVolm;
    unsigned int LocalEntitySizes[4] = {0,0,0,0};
    unsigned int LocalEntityStarts[4] = {0,0,0,0};
    int nprocs, procID;

  public:

    template<typename GlobalIterator, typename LocalIterator>
    void AddToMapIteratorIfUnique(std::map<GlobalIterator,LocalIterator> EntityMap, GlobalIterator I, LocalIterator J){
      if( EntityMap.find(I) == EntityMap.end() ){
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
          if( node.processor_id() == LProcID){
            if( elem->is_vertex(nodeID)   ) AddToMapIteratorIfUnique<unsigned int, unsigned int>(Global_to_LVert, nodeID, LocalEntitySizes[0]);
            if( elem->is_edge(nodeID)     ) AddToMapIteratorIfUnique<unsigned int, unsigned int>(Global_to_LEdge, nodeID, LocalEntitySizes[1]);
            if( elem->is_face(nodeID)     ) AddToMapIteratorIfUnique<unsigned int, unsigned int>(Global_to_LFace, nodeID, LocalEntitySizes[2]);
            if( elem->is_internal(nodeID) ) AddToMapIteratorIfUnique<unsigned int, unsigned int>(Global_to_LVolm, nodeID, LocalEntitySizes[3]);
          }
        }
      }

      //Find the entity sizes on each processor
      //and 
      std::vector<unsigned int> procEntitySizesGlobal;
      procEntitySizesGlobal.clear();
      for(int I=0; <4*nprocs; I++) procEntitySizesGlobal.push_back(0);


      procEntitySizesGlobal[procID*4 + 0] = LocalEntitySizes[0];
      procEntitySizesGlobal[procID*4 + 1] = LocalEntitySizes[1];
      procEntitySizesGlobal[procID*4 + 2] = LocalEntitySizes[2];
      procEntitySizesGlobal[procID*4 + 3] = LocalEntitySizes[3];
      MPI_Allreduce(&procEntitySizesGlobal.front(), &procEntitySizesGlobal.front(), &procEntitySizesGlobal.size()
                  , MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

      if(procID != 0){
        for(int I=0; I<procID; I++){
          int K = ;
          LocalEntityStarts[0] += procEntitySizesGlobal[I*4 + 0];
          LocalEntityStarts[1] += procEntitySizesGlobal[I*4 + 1];
          LocalEntityStarts[2] += procEntitySizesGlobal[I*4 + 2];
          LocalEntityStarts[3] += procEntitySizesGlobal[I*4 + 3];
        }
      }
    };

    bool Is_LocalEdge(unsigned int nodeID){
      if( Global_to_LEdge.find(nodeID) == Global_to_LEdge.end() ) return false;
      return false;
    }
	
    int EdgeLocalID(unsigned int nodeID){
      return Global_to_LEdge[nodeID];
    }
	//Nothing Interesting
};
