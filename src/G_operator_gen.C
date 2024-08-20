#include "G_operator_gen.hpp"

//The class constructor
G_operator::G_operator(EquationSystems & es, SupplementaryEntityIDs & SupEiDs){
  if(is_parallel){
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &procID);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  }
  Make_Edge_Map(es);
  Set_G_Operator();
};

//Cantor counting function
unsigned int G_operator::Cantors_Counter(unsigned int I , unsigned int J){
  return (((I+J)*(I+J+1))/2) + J;
};

std::pair<unsigned int, unsigned int> G_operator::Cantors_CounterInv(unsigned int K){
  float K_float = float(K);
  unsigned int W = std::floor(0.5*(std::sqrt(8.0*K_float + 1.0) - 1.0));
  unsigned int T = (W*W + W)/2;
  return make_pair( (K-T), (W-K+T) );
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
      std::pair<unsigned int, unsigned int> edge_key;
      if(m < n) edge_key = make_pair(m,n);
      if(n < m) edge_key = make_pair(n,m);
      edge_map[edge_key] = Cantors_Counter(edge_key.first,edge_key.second) ;
    }
  }


  ntot_edges_local = edge_map.size();  
}


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