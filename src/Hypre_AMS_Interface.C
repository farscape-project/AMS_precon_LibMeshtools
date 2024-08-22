#include "Hypre_AMS_Interface.hpp"

//The class constructor
Hypre_AMS_Interface::Hypre_AMS_Interface(EquationSystems & es): _SupEiDs(es){
  if(is_parallel){
    ierr = MPI_Comm_rank(mesh.comm(), &procID);
    ierr = MPI_Comm_size(mesh.comm(), &nprocs);
  }
  Make_Edge_Map(es);
  Set_Hypre_AMS_Interface();
};


// Makes the edge map
void Hypre_AMS_Interface::Make_Edge_Map(EquationSystems & es)
{
  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();


  // Form the edge-node pairing
  for (const auto & elem : mesh.active_local_element_ptr_range())
  {
    const unsigned int nedges = elem->n_edges();
    for(unsigned int I=0; I<nedges; I++)
    {
      //Find global nodeIDs of the edge endpoints
      unsigned int EdgeID = elem->node_id(elem->edge_nodes_map[I][2]))
      if( _SupEiDs->Is_LocalEdge(EdgeID) ){
        //Form a pair of the to vertices at either end of the edge
		unsigned int EdgeLocalID = EdgeLocalID(EdgeID);
		unsigned int m = elem->node_id(elem->edge_nodes_map[I][0]));
        unsigned int n = elem->node_id(elem->edge_nodes_map[I][1]));
        std::pair<unsigned int, unsigned int> edge_EndPair;
        edge_EndPair = make_pair(m,n);
        edge_map[EdgeLocalID] = edge_EndPair;
      }
    }
  }
  ntot_edges_local = edge_map.size();  
}


// Sets the G-operator matrix using the PETSc-hypre 
// interface using the IJ matrix interface
void Hypre_AMS_Interface::Set_Hypre_AMS_Interface(EquationSystems & es){

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  //size up and set up the arrays
  int nrows;
  int *ncols, *rows, *cols;
  double *Matvalues, *Vecvalues;

  ilower = _SupEiDs->LocalEntityStarts[1];         //local lower bound for global edge number
  iupper = ilower + _SupEiDs->LocalEntitySizes[1]; //local upper bound for global edge number
  jlower = ;                                       //local lower bound for global vertex number
  jupper = ;                                       //local lower bound for global vertex number


  //Set the sizing aray values
  nrows  =  ProcEdgeSize[procID]
  ncols  = new int[nrows];
  rows   = new int[nrows];
  for(int I=0; I<nrows; I++){
    ncols[I] = 2;
    row[I] = I + ilower;
  }


  // Set the matrix values and the 
  // from the edge-map
  //(There are only two column entries per row
  // so no advanced calculations are really
  // needed for this)
  cols   = new int[2*nrows]; 
  Matvalues = new double[2*nrows];


  // Iterator for the edge-map
  int K=0;
  std::map<int,std::pair<unsigned int, unsigned int>>::iterator it;
  for(it = edge_map.begin(); it != edge_map.end(); it++){
    //Assign to CSR matrix+value
    cols[K] = it.first;
	Matvalues[K] =  1.0;
    K++;
    cols[K] = it.second;
	Matvalues[K] = -1.0;
    K++;
  };

  //coordinates at vertices (PETSc Vector)
  Vec  par_xcoord, par_ycoord, par_zcoord; 
  petscErr = VecCreate(mesh.comm(),&par_xcoord);
  petscErr = VecCreate(mesh.comm(),&par_ycoord);
  petscErr = VecCreate(mesh.comm(),&par_zcoord)

  //Set the coordinate vector sizes and paritions
  int CoordsSize = _SupEiDs->Global_to_LVert.size();
  petscErr = VecSetSizes(par_xcoord,PETSC_DECIDE,CoordsSize);
  petscErr = VecSetFromOptions(par_xcoord);
  petscErr = VecDuplicate(par_xcoord,&par_ycoord);
  petscErr = VecDuplicate(par_xcoord,&par_zcoord);

  //Setting the vector-coordinate Values
  int istart,iend;
  VecGetOwnershipRange(par_xcoord,&istart,&iend);
  std::map<unsigned int, unsigned int>::iterator it;
  it = _SupEiDs->Global_to_LVert.begin();
  for(int I=istart; I<iend; I++){
    int nodeID = it->first;
    PetscScalar x = (PetscScalar)( mesh. );
    PetscScalar y = (PetscScalar)( mesh. );
    PetscScalar z = (PetscScalar)( mesh. );
    VecSetValues(par_xcoord,1,&I,&x,INSERT_VALUES);
    VecSetValues(par_ycoord,1,&I,&y,INSERT_VALUES);
    VecSetValues(par_zcoord,1,&I,&z,INSERT_VALUES);
	it++;
  }


  //Create the empty matrix and vectors
  petscErr = MatCreate(mesh.comm(), par_G);


  //Set the G-Operator matrix
  petscErr = PCHYPRESetDiscreteGradient(PC pc, par_G);


  //Multiply the G-operator by the coordinates
  //to get the edge vectors
  petscErr = MatMult(par_G, par_xcoord, par_xvec);
  petscErr = MatMult(par_G, par_ycoord, par_yvec);
  petscErr = MatMult(par_G, par_zcoord, par_zvec);


  //Set the G-operator matrix
  petscErr = PCHYPRESetEdgeConstantVectors(PC pc, par_xvec, par_yvec, par_zvec);


  //Clean-up the extra arrays
  ierr = VecDestroy(&par_xcoord);
  ierr = VecDestroy(&par_ycoord);
  ierr = VecDestroy(&par_zcoord);
  delete[] ncols, rows, cols, values;
};