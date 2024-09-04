#include "Hypre_AMS_Interface.h"

//The class constructor
Hypre_AMS_Interface::Hypre_AMS_Interface(EquationSystems & es){
  _SupEiDs.FormEntityMaps(es);
  _SupEiDs.FormVertexMaps(es);
  Allocate_G_Operator(es);
  Make_Edge_Map(es);
  Set_Hypre_AMS_Interface(es);
};

void Hypre_AMS_Interface::Allocate_G_Operator(EquationSystems & es)
{
  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  //Create the empty matrix and vectors
  petscErr = MatCreate(mesh.comm().get(), &par_G);
  petscErr = MatSetType(par_G, MATMPIAIJ); 
  petscErr = MatSetSizes(par_G, _SupEiDs.local_num_rows, _SupEiDs.local_num_cols, _SupEiDs.total_num_rows, _SupEiDs.total_num_cols); 
  PetscInt d_nz = 2;
  PetscInt o_nz = 2;
  petscErr = MatMPIAIJSetPreallocation(par_G, d_nz, NULL, o_nz, NULL);
}


// Makes the edge map
void Hypre_AMS_Interface::Make_Edge_Map(EquationSystems & es)
{
  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  std::vector<PetscScalar> vals{1, -1};
  unsigned int first_index, second_index;
  unsigned int first_col_index, second_col_index;

  // Form the edge-node pairing
  for (const auto & elem : mesh.active_local_element_ptr_range())
  {
    const unsigned int nedges = elem->n_edges();
    for(unsigned int I=0; I<nedges; I++)
    {
      //Find global nodeIDs of the edge endpoints
      unsigned int EdgeID = elem->node_id(elem->local_edge_node(I, 2));
      const Node & EdgeRef = mesh.node_ref(EdgeID);
      if( _SupEiDs.Is_LocalEdge(EdgeID) ){
        //Form a pair of the to vertices at either end of the edge
		    unsigned int EdgeLocalID = _SupEiDs.EdgeLocalID(EdgeID);
		    unsigned int m = elem->node_id(elem->local_edge_node(I, 0));
        unsigned int n = elem->node_id(elem->local_edge_node(I, 1));
        const Node & node_m = mesh.node_ref(m);
        const Node & node_n = mesh.node_ref(n);
        const short sign = node_m > node_n ? 1 : -1;
        if(sign > 0) {first_index = m; second_index = n;} else {first_index = n; second_index = m;}
        //std::vector<PetscInt> row_index{EdgeLocalID+_SupEiDs.LocalEntityStarts[1]};
        std::vector<PetscInt> row_index{EdgeRef.dof_number(0,0,0)};
        first_col_index = _SupEiDs.GlobalVertexID_to_SeqID[first_index];
        second_col_index = _SupEiDs.GlobalVertexID_to_SeqID[second_index];
        std::vector<PetscInt> col_index{first_col_index , second_col_index};
        MatSetValues(par_G,  1, row_index.data(),  2, col_index.data(), vals.data(), INSERT_VALUES);
      }
    }
  }

  petscErr = MatAssemblyBegin(par_G, MAT_FINAL_ASSEMBLY);
  petscErr = MatAssemblyEnd(par_G, MAT_FINAL_ASSEMBLY); 
 
  //MatView(par_G, PETSC_VIEWER_STDOUT_WORLD);
}


// Sets the G-operator matrix using the PETSc-hypre 
// interface using the IJ matrix interface
void Hypre_AMS_Interface::Set_Hypre_AMS_Interface(EquationSystems & es){

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  //coordinates at vertices (PETSc Vector)
  Vec  par_xcoord, par_ycoord, par_zcoord; 
  petscErr = VecCreateMPI(mesh.comm().get(),_SupEiDs.local_num_cols,_SupEiDs.total_num_cols,&par_xcoord);
  petscErr = VecCreateMPI(mesh.comm().get(),_SupEiDs.local_num_cols,_SupEiDs.total_num_cols,&par_ycoord);
  petscErr = VecCreateMPI(mesh.comm().get(),_SupEiDs.local_num_cols,_SupEiDs.total_num_cols,&par_zcoord);
 
  for(auto it = _SupEiDs.Global_to_LVert.begin(); it != _SupEiDs.Global_to_LVert.end(); it++){
    int nodeID = it->first;
    const Node & node = mesh.node_ref(nodeID);
    PetscScalar x = (PetscScalar)( node(0) );
    PetscScalar y = (PetscScalar)( node(1) );
    PetscScalar z = (PetscScalar)( node(2) );
    int localID = it->second;
    PetscInt I = _SupEiDs.LocalEntityStarts[0]+localID;
    VecSetValues(par_xcoord,1,&I,&x,INSERT_VALUES);
    VecSetValues(par_ycoord,1,&I,&y,INSERT_VALUES);
    VecSetValues(par_zcoord,1,&I,&z,INSERT_VALUES);
  }

  petscErr = VecAssemblyBegin(par_xcoord); 
  petscErr = VecAssemblyEnd(par_xcoord); 

  petscErr = VecAssemblyBegin(par_ycoord); 
  petscErr = VecAssemblyEnd(par_ycoord); 

  petscErr = VecAssemblyBegin(par_zcoord); 
  petscErr = VecAssemblyEnd(par_zcoord); 

  //VecView(par_xcoord, PETSC_VIEWER_STDOUT_WORLD);


/*

  //Set the G-Operator matrix
  //petscErr = PCHYPRESetDiscreteGradient(pc, par_G);


  //Multiply the G-operator by the coordinates
  //to get the edge vectors
  petscErr = MatMult(par_G, par_xcoord, par_xvec);
  petscErr = MatMult(par_G, par_ycoord, par_yvec);
  petscErr = MatMult(par_G, par_zcoord, par_zvec);


  //Set the G-operator matrix
  //petscErr = PCHYPRESetEdgeConstantVectors(pc, par_xvec, par_yvec, par_zvec);


  //Clean-up the extra arrays
  int ierr = VecDestroy(&par_xcoord);
  ierr = VecDestroy(&par_ycoord);
  ierr = VecDestroy(&par_zcoord);
  //delete[] ncols, rows, cols, Matvalues, Vecvalues;
  */
};