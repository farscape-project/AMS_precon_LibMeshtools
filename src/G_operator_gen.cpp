#include "G_operator_gen.hpp"

//The class constructor
G_operator::G_operator(EquationSystems & es, const std::string & system_name){
  Make_Edge_Map(es, system_name);
  Size_G_Operator();
  Set_G_Operator();
};


// Makes the edge map (these are non-unique)
void G_operator::Make_Edge_Map(EquationSystems & es, const std::string & system_name)
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
}


// Sizes up the G-operator matrix for the PETSc-hypre 
// interface
void G_operator::Size_G_Operator(){
  //=====
  // Get the MPI jobsize
  //=====
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &procID);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


//===================================================================
// This is an extremely wasteful procedure (stores and nproc size 
// array on every process) however i am using it as realistically 
// I am not dealing with exascale (yet) an even then an int array
// is not too bad
//===================================================================
  //=====
  // Set the edge size
  //=====
  ProcEdgeSize.clear();
  for(int I=0; I<nprocs; I++) ProcEdgeSize.push_back(0);

  //=====
  // Find the global number of edges 
  //=====
  ProcEdgeSize[procID] = ntot_edges_local;
  MPI_Allreduce(&ProcEdgeSize.front(), &ProcEdgeSize.front(), &ProcEdgeSize.size(), MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  ntot_edges_global = 0; //Just making sure its zero;
  for(int I=0; I<nprocs; i++) ntot_edges_global = ntot_edges_global + ProcEdgeSize[I];

  //=====
  // Find the ilower and iupper
  // (local Edge bounds)
  //=====
  ilower = 0;
  for(int I=0; I<procID; I++) ilower = ilower + ProcEdgeSize[I];
  iupper = ilower + ProcEdgeSize[procID];
//===================================================================
//===================================================================

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
  for(it = edge_map.begin(); it != edge_map.end(); it++){
    cols[K] = it.first;
	values[K] =  1.0;
    K++;
    cols[K] = it.second;
	values[K] = -1.0;
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