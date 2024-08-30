#include "ElementNode_EntityMapper.h"

//Form the entity Maps
void SupplementaryEntityIDs::FormEntityMaps(EquationSystems & es){
  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  //std::cout << "MyProcID " << procID  << " We are HERE in SupplementaryEntityIDs::FormEntityMaps" << std::endl; 

  //Find maps entities in global node numberings from local process
  //to a unique contiguous ID 
  unsigned int LProcID = mesh.processor_id();

  std::cout << " LProcID " << mesh.processor_id() << " MyProcID " << procID << std::endl;

  for (const auto & elem : mesh.active_local_element_ptr_range()){
    for(unsigned int I=0; I<elem->n_nodes(); I++){
      unsigned int nodeID = elem->node_id(I);

      const Node & node = mesh.node_ref(nodeID);
      if( node.processor_id() == LProcID){
            
        if( elem->is_vertex(I)   ) 
        {
       // std::cout << "nodeID " << nodeID << std::endl;
        AddToMapIteratorIfUnique<unsigned int, unsigned int>(Global_to_LVert, nodeID, LocalEntitySizes[0], elem, "NODE");
       // std::cout << "LOOP LocalEntitySizes[0] " << LocalEntitySizes[0] << std::endl;
        }
        if( elem->is_edge(I)     ) AddToMapIteratorIfUnique<unsigned int, unsigned int>(Global_to_LEdge, nodeID, LocalEntitySizes[1], elem, "EGDE");
        if( elem->is_face(I)     ) AddToMapIteratorIfUnique<unsigned int, unsigned int>(Global_to_LFace, nodeID, LocalEntitySizes[2], elem, "FACE");
        if( elem->is_internal(I) ) AddToMapIteratorIfUnique<unsigned int, unsigned int>(Global_to_LVolm, nodeID, LocalEntitySizes[3], elem, "INTERNAL");
      }
    }
  }

  //Find the entity sizes on each processor
  //and 
  // std::cout << "ONE LocalEntitySizes[0] " << LocalEntitySizes[0] << std::endl;

  std::vector<unsigned int> procEntitySizesGlobal;
  procEntitySizesGlobal.clear();
  for(int I=0; I<4*nprocs; I++) procEntitySizesGlobal.push_back(0);

  // std::cout << "TWO LocalEntitySizes[0] " << LocalEntitySizes[0] << std::endl;


  procEntitySizesGlobal[procID*4 + 0] = LocalEntitySizes[0];
  procEntitySizesGlobal[procID*4 + 1] = LocalEntitySizes[1];
  procEntitySizesGlobal[procID*4 + 2] = LocalEntitySizes[2];
  procEntitySizesGlobal[procID*4 + 3] = LocalEntitySizes[3];
  MPI_Allreduce(&procEntitySizesGlobal.front(), &procEntitySizesGlobal.front(), procEntitySizesGlobal.size()
              , MPI_UNSIGNED, MPI_SUM, mesh.comm().get());

   //std::cout << "THREE LocalEntitySizes[0] " << LocalEntitySizes[0] << std::endl;

  if(procID != 0){
    for(int I=0; I<procID; I++){
      LocalEntityStarts[0] += procEntitySizesGlobal[I*4 + 0];
      LocalEntityStarts[1] += procEntitySizesGlobal[I*4 + 1];
      LocalEntityStarts[2] += procEntitySizesGlobal[I*4 + 2];
      LocalEntityStarts[3] += procEntitySizesGlobal[I*4 + 3];
      //std::cout << "procID " << procID << " THREE LocalEntitySizes[0] " << LocalEntitySizes[0] << std::endl;
    }
  }

  std::cout << "MyProcID " << procID  << " LocalEntitySizes[0] "  << LocalEntitySizes[0] << std::endl; 
  std::cout << "MyProcID " << procID  << " LocalEntitySizes[1] " << LocalEntitySizes[1] << std::endl; 
  std::cout << "MyProcID " << procID  << " LocalEntitySizes[2] "  << LocalEntitySizes[2] << std::endl; 
  std::cout << "MyProcID " << procID  << " LocalEntitySizes[3] " << LocalEntitySizes[3] << std::endl; 
};


bool SupplementaryEntityIDs::Is_LocalEdge(unsigned int nodeID){
  if( Global_to_LEdge.find(nodeID) == Global_to_LEdge.end() ) return false;
  return false;
}
	
int SupplementaryEntityIDs::EdgeLocalID(unsigned int nodeID){
  if( Global_to_LEdge.find(nodeID) == Global_to_LEdge.end() ) return Global_to_LEdge[nodeID];
  return -1;
}

