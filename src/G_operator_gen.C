#include "G_operator_gen.hpp"

//The class constructor
G_operator::G_operator(EquationSystems & es){
  if(is_parallel){
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &procID);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  }
  Make_Edge_Map(es);
  Size_G_Operator();
  Set_G_Operator();
};

//Cantor counting function
int G_operator::Cantors_Counter(unsigned int I , unsigned int J){
  return (((I+J)*(I+J+1))/2) + J;
};

std::pair<unsigned int, unsigned int> G_operator::Cantors_CounterInv(int K){


  return;
};


void G_operator::Map_Proc_Neighbors(EquationSystems & es){
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
      edge_map[edge_key] = Cantors_Counter(edge_key.first,edge_key.second) ;
    }
  }
}


// Makes the edge map
void G_operator::Make_Edge_Map(EquationSystems & es)
{
  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  //=====
  // Form the initial edge-node-pair-map
  //=====
  // For serial applications this is sufficient
  // for calculating G-operator matrix needed in
  // AMS however in parallel the additional step of
  // removing remote duplicates is needed.
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
      edge_map[edge_key] = Cantors_Counter(edge_key.first,edge_key.second) ;
    }
  }

  //Remove the non-local duplicates
  if(is_parallel){
    prune_Remote_Duplicate_Edges();
  }
  ntot_edges_local = edge_map.size();  
}

// Called within make edgemap if the system is parallel
// removes all remote edges and assigns unique edges
// a unique proc ID, This is done by large neighborhood 
// search
void G_operator::prune_Remote_Duplicate_Edges(){

  //Forming an extended neighborhood to search
  std::map<unsigned int, unsigned int> EdgeOwnerMap;
  std::map<unsigned int, unsigned int>::iterator jt;
  std::vector<unsigned int> LECantorIters, TLECantorIters, IterStart; 


  //Find out their edge sizes and allocate the space on the vectors
  int K=0;
  IterStart.clear();
  TLECantorIters.clear();
  for(int I=0; I<ProcEdgeSize.size(); I++){
    IterStart.push_back(K);
    for(int J=0; J<ProcEdgeSize[I]; J++){
      TLECantorIters.push_back(0);
      K++;
    }
  }


  // Send and recieve the Cantor ID's vector 
  // (for the edges) to neighbouring processors
  MPI_Request requests[nmessages];
  MPI_Status  statuses[nmessages];
  for(int I=0; I<nmessages; I++) //sends local edge-cantor-iterators to neighbors
    ier = MPI_Isend(&LECantorIters.front(),LECantorIters.size(),MPI_UNSIGNED, ProcNeighbors[I], procID
                    ,MPI_COMM_WORLD, &requests[I]);

  for(int I=0; I<nmessages; I++) //recieves non-local edge-cantor-iterators to neighbors
    ier = MPI_Recv(&TLECantorIters[IterStart[I]], ProcEdgeSize[I], MPI_UNSIGNED, ProcNeighbors[I], ProcNeighbors[I]
                 , MPI_COMM_WORLD, &statuses[I]);

  //Just waiting for all communications at this phase to finish
  ier = MPI_Waitall(nmessages, requests, statuses)


  // Form a map of the processor owner of each of the edges
  // the heuristic used here is to make the processor containing
  // the edge with the minimum procID the owner of the edge
  for(int I=0; I<ProcEdgeSize.size(); I++){
    for(int J=IterStart[I]; J<ProcEdgeSize[I]; J++){
      int K = 0; // The cantor iterator of edge pair
      int L = 0; // The processor ID of that entry
      if(edge_map.find(edgeKey) ==  edge_map.end() ){
        EdgeOwnerMap[edgeKey] = L;
      }else{
        EdgeOwnerMap[edgeKey] = std::min(EdgeOwnerMap[edgeKey],L);
      }
    }
  }

  //Remove any local edges that do not/no-longer belong to the
  //local process
  for(jt = EdgeOwnerMap.begin(); jt != EdgeOwnerMap.end(); jt++){
    std::pair<unsigned int, unsigned int>  edgeKey = Cantors_CounterInv(jt.first);
    if(edge_map.find(edgeKey) !=  edge_map.end() ){ //Checks if edge is local
      if(jt.second != procID){ //Checks if edge belongs to local process
        edge_map.erase(edgeKey); //Erase non-locally owned edges
      }
    }
  }
  //Successfully pruned all non-local edges 
  //now just cleanup all the search arrays
  LECantorIters.clear();
  TLECantorIters.clear();
  IterStart.clear();
};


// Sizes up the G-operator matrix for the PETSc-hypre 
// interface
void G_operator::Size_G_Operator(){

  //=====
  // Set the edge size
  //=====
  ProcEdgeSize.clear();
  for(int I=0; I<nprocs; I++) ProcEdgeSize.push_back(0);

  //=====
  // Find the global number of edges 
  //=====
  ProcEdgeSize[procID] = ntot_edges_local;
  MPI_Allreduce(&ProcEdgeSize.front(), &ProcEdgeSize.front(), &ProcEdgeSize.size()
              , MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

  ntot_edges_global = 0; //Just making sure its zero;
  for(int I=0; I<nprocs; i++) ntot_edges_global = ntot_edges_global + ProcEdgeSize[I];

  //=====
  // Find the ilower and iupper
  // (local Edge bounds)
  //=====
  ilower = 0;
  for(int I=0; I<procID; I++) ilower = ilower + ProcEdgeSize[I];
  iupper = ilower + ProcEdgeSize[procID];


  //=====
  // Find the minimum and maximum global
  // node numbers of local process
  // (This is a wasteful search however
  // its only done once)
  //=====
  //declare and initialise the nodeIDS
  unsigned int minLNodeID, maxLNodeID;
  maxLNodeID = 0;
  minLNodeID = edge_map[0].first;
  for(it = edge_map.begin(); it != edge_map.end(); it++){
    //Check nodes in edges and find max
    maxLNodeID = std::max(maxLNodeID,it.first);
    maxLNodeID = std::max(maxLNodeID,it.second);

    //Check nodes in edges and find min
    minLNodeID = std::min(minLNodeID,it.first);
    minLNodeID = std::min(minLNodeID,it.second);
  }

  //=====
  // Find the jlower and jupper
  // (local node bounds)
  //=====
  jlower = minLNodeID;
  jupper = maxLNodeID;
};


// Sets the G-operator matrix using the PETSc-hypre 
// interface using the IJ matrix interface
void G_operator::Set_G_Operator(){
  int nrows;
  int *ncols, *rows, *cols;
  double *values;

  //=====
  //Set the sizing aray values
  //=====
  nrows  =  ProcEdgeSize[procID]
  ncols  = new int[nrows];
  rows   = new int[nrows];
  for(int I=0; I<nrows; I++){
    ncols[I] = 2;
    row[I] = I + ilower;
  }

  //=====
  //Set the matrix values and the 
  // from the edge-map
  //=====
  //(There are only two entries per row
  //so no advanced calculations are really
  //needed for this)
  cols   = new int[2*nrows]; 
  values = new double[2*nrows];
  int K=0;

  // Iterator for the edge-map
  std::map<int,std::pair<unsigned int, unsigned int>>::iterator it;
  for(it = edge_map.begin(); it != edge_map.end(); it++){
    //Find orientation of the edge
    double signAB=1.0;
    signAB = signAB*InnerProductSign<double>(,)
    signAB = signAB*InnerProductSign<double>(,)
    signAB = signAB*InnerProductSign<double>(,)

    //Assign to CSR matrix+value
    cols[K] = it.first;
	values[K] =  1.0*signAB;
    K++;
    cols[K] = it.second;
	values[K] = -1.0*signAB;
    K++;
  };

  //=====
  //Generate the matrix
  //=====
  HYPRE_IJMatrixCreate(comm, ilower, iupper, jlower, jupper, &par_G_ij);
  HYPRE_IJMatrixSetObjectType(par_G_ij, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(par_G_ij);

  //=====
  //Set matrix coefficients
  //=====
  HYPRE_IJMatrixSetValues(par_G_ij, nrows, ncols, rows, cols, values);
  HYPRE_IJMatrixAssemble(par_G_ij);
  HYPRE_IJMatrixGetObject(par_G_ij, (void **) &par_G);

  //=====
  //Clean-up the extra arrays
  //=====
  delete[] ncols, rows, cols, values;
};