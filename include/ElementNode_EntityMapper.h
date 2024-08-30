#ifndef ELEMENTNODE_ENTITYMAPPER_H
#define ELEMENTNODE_ENTITYMAPPER_H 

#include <map>
#include "libmesh/equation_systems.h"
#include "libmesh/elem.h"

using namespace libMesh;

class SupplementaryEntityIDs
{
  public:
    //These maps store a global contiguous numbering of element subentities
    //And map these to the global Node/Dof numbering system within LibMesh
    std::map<unsigned int, unsigned int> Global_to_LVert;
    std::map<unsigned int, unsigned int> Global_to_LEdge;
    std::map<unsigned int, unsigned int> Global_to_LFace;
    std::map<unsigned int, unsigned int> Global_to_LVolm;
    unsigned int LocalEntitySizes[4] = {0,0,0,0};
    unsigned int LocalEntityStarts[4] = {0,0,0,0};
    processor_id_type nprocs{global_n_processors()};
    processor_id_type procID{global_processor_id()};


    //Template that increases the size of a map
	//if the new entry is unique
    template<typename GlobalIterator, typename LocalIterator>
    void AddToMapIteratorIfUnique(std::map<GlobalIterator,LocalIterator> EntityMap, GlobalIterator I, LocalIterator J){
      if( EntityMap.find(I) == EntityMap.end() ){
        EntityMap[I] = J;
        J++;
      }
    }

    //Form the entity Maps
    void FormEntityMaps(EquationSystems & es);

    //Checks whether the global nodeID provided
	//corresponds to an edge that is owned by the
	//local process
    bool Is_LocalEdge(unsigned int nodeID);

    //Returns the local contiguous EdgeID 
	//when given a valid global nodeID
	//Otherwise returns a value of -1
    int EdgeLocalID(unsigned int nodeID);
};
#endif