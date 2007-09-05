#include "Output.h"

Output::Output(Grid *g) {

  grid = g;

}

void Output::printGMSGrid() {

  FILE *fp = NULL;
  
  PetscPrintf(PETSC_COMM_WORLD,"Printing GMS grid.\n");
  PetscFOpen(PETSC_COMM_WORLD,"gms.3dg","w",&fp);
  PetscFPrintf(PETSC_COMM_WORLD,fp,"GRID3D\n");
  PetscFPrintf(PETSC_COMM_WORLD,fp,"TYPE 1\n");
  PetscFPrintf(PETSC_COMM_WORLD,fp,"IJK +y +x +z\n");
  PetscFPrintf(PETSC_COMM_WORLD,fp,"NUMBERING 0\n");
  double *ptr = grid->getOriginPtr();
  PetscFPrintf(PETSC_COMM_WORLD,fp,"ORIGIN %.8e %.8e %.8e\n",ptr[0],ptr[1],ptr[2]);
  ptr = NULL;
  PetscFPrintf(PETSC_COMM_WORLD,fp,"ROTZ %f\n",grid->getRotationDegrees());
  PetscFPrintf(PETSC_COMM_WORLD,fp,"DIM %d %d %d\n",grid->getNx()+1,
               grid->getNy()+1,grid->getNz()+1);

  double sum = 0.;
  PetscFPrintf(PETSC_COMM_WORLD,fp,"%.8e\n",sum);
  for (int i=0; i<grid->getNx(); i++) {
    sum += grid->getDx(i);
    PetscFPrintf(PETSC_COMM_WORLD,fp,"%.8e\n",sum);
  }
  sum = 0.;
  PetscFPrintf(PETSC_COMM_WORLD,fp,"%.8e\n",sum);
  for (int j=0; j<grid->getNy(); j++) {
    sum += grid->getDy(j);
    PetscFPrintf(PETSC_COMM_WORLD,fp,"%.8e\n",sum);
  }
  sum = 0.;
  PetscFPrintf(PETSC_COMM_WORLD,fp,"%.8e\n",sum);
  for (int k=0; k<grid->getNz(); k++) {
    sum += grid->getDz(k);
    PetscFPrintf(PETSC_COMM_WORLD,fp,"%.8e\n",sum);
  }
  PetscFPrintf(PETSC_COMM_WORLD,fp,"MAT\n");
  Vec v = grid->getGridCellMaterialIDs();
  writeIntVectorInNaturalOrder(fp,v,0);
  PetscFPrintf(PETSC_COMM_WORLD,fp,"ACTIVE\n");
  v = grid->getGridCellActivities();
  writeIntVectorInNaturalOrder(fp,v,0);
  VecDestroy(v);

  PetscFClose(PETSC_COMM_WORLD,fp);

}

void Output::printGMSDataSet(char *filename, Vec v) {
    
  FILE *fp = NULL;
  
  PetscPrintf(PETSC_COMM_WORLD,"Printing GMS data set (%s).\n",filename);
  PetscFOpen(PETSC_COMM_WORLD,filename,"w",&fp);
  PetscFPrintf(PETSC_COMM_WORLD,fp,"DATASET\n");
  PetscFPrintf(PETSC_COMM_WORLD,fp,"OBJTYPE \"grid3d\"\n");
  PetscFPrintf(PETSC_COMM_WORLD,fp,"BEGSCL\n");
  PetscFPrintf(PETSC_COMM_WORLD,fp,"ND %d\n",grid->getN());
  PetscFPrintf(PETSC_COMM_WORLD,fp,"NC %d\n",grid->getN());
  char word[32];
  int i;
  for (i=0; i<32; i++) {
    if (filename[i] == 46) break; // 46 = "."
    word[i] = filename[i];
  }
  word[i] = '\0';
  PetscFPrintf(PETSC_COMM_WORLD,fp,"NAME %s\n",word);
  PetscFPrintf(PETSC_COMM_WORLD,fp,"TS 0 0.\n");
  writeIntVectorInNaturalOrder(fp,v,1);
  PetscFPrintf(PETSC_COMM_WORLD,fp,"ENDDS\n");

  PetscFClose(PETSC_COMM_WORLD,fp);
}

void Output::writeIntVectorInNaturalOrder(FILE *fp, Vec v, int one_per_line) {
  Vec natural;
  PetscScalar *v_ptr = NULL;
  grid->getVectorNatural(&natural);
  grid->globalToNatural(v,natural);
  int max_cells;
  MPI_Allreduce(&grid->num_cells_local,&max_cells,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);
  int *values = new int[max_cells];
  VecGetArray(natural,&v_ptr);
  for (int i=0; i<grid->num_cells_local; i++) 
    values[i] = int(v_ptr[i]);
  VecRestoreArray(natural,&v_ptr);
  VecDestroy(natural);
  if (myrank == 0) {
    int count = 0;
    for (int i=0; i<grid->num_cells_local; i++) {
      if (one_per_line) {
        PetscFPrintf(PETSC_COMM_WORLD,fp,"%d\n",values[i]);
      }
      else {
        PetscFPrintf(PETSC_COMM_WORLD,fp,"%d ",values[i]);
        if (++count%20 == 0) PetscFPrintf(PETSC_COMM_WORLD,fp,"\n");
      }
    }
    for (int iproc=1; iproc<commsize; iproc++) {
      MPI_Status status;
      MPI_Probe(iproc,MPI_ANY_TAG,PETSC_COMM_WORLD,&status);
      int nrecv = status.MPI_TAG;
      MPI_Recv(values,nrecv,MPI_INTEGER,iproc,MPI_ANY_TAG,PETSC_COMM_WORLD,&status);
      for (int i=0; i<nrecv; i++) {
        if (one_per_line) {
          PetscFPrintf(PETSC_COMM_WORLD,fp,"%d\n",values[i]);
        }
        else {
          PetscFPrintf(PETSC_COMM_WORLD,fp,"%d ",values[i]);
          if (++count%20 == 0) PetscFPrintf(PETSC_COMM_WORLD,fp,"\n");
        }
      }
    }
    if (!one_per_line && count%10 != 0) PetscFPrintf(PETSC_COMM_WORLD,fp,"\n");
  }
  else {
    MPI_Send(values,grid->num_cells_local,MPI_INTEGER,0,grid->num_cells_local,
             PETSC_COMM_WORLD);
  }
  delete(values);
}

void Output::printBoundarySets() {
  Vec v;
  PetscScalar *v_ptr = NULL;
  grid->getVectorGlobal(&v);
  BoundarySet *cur_set = grid->boundary_sets;  
  while (cur_set) {
    VecGetArray(v,&v_ptr);
    for (int i=0; i<grid->num_cells_local; i++)
      v_ptr[i] = 0.;
    Connection *cur_connection = cur_set->list;
    while (cur_connection) {
      v_ptr[cur_connection->cell] += 1;
      cur_connection = cur_connection->next;
    }
    VecRestoreArray(v,&v_ptr);
    char filename[32];
    strcpy(filename,cur_set->name);
    strcat(filename,".dat");
    printGMSDataSet(filename,v);
    cur_set = cur_set->next;
  }
  VecDestroy(v);
}

void Output::printHDFMesh() {
#if 0
  char filename[32];
  strcpy(filename,"grid.h5");

  hid_t file_id; // file identifier
  hid_t prop_id; // property list
  hid_t cell_grp_id; // cell group identifier
  hid_t vert_grp_id; // vertex group identifier
  hid_t data_space_id; // data space identifier (file space)
  hid_t data_set_id; // data set identifier (memory space)
  herr_t status;

  // set up parallel access, if applicable
  prop_id = H5Pcreate(H5P_FILE_ACCESS);
#ifdef PARALLEL
  H5Pset_fapl_mpio(prop_id,PETSC_COMM_WORLD,MPI_INFO_NULL);
#endif

// Create a new file collectively and release property list identifier
  /* Access flags:
     H5F_ACC_TRUNC  - overwrites existing file, deleting original contents
     H5F_ACC_EXCL   - open should fail if file exist
     H5F_ACC_RDONLY - read access only
     H5F_ACC_RDWR   - read and write access
  */                                          // create id // access id
  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, prop_id);
  H5Pclose(prop_id);

// Open cell group  
/*
  herr_t (*old_func)(void*);
  void  *old_client_data;
  H5Eget_auto(&old_func, &old_client_data);
  H5Eset_auto(NULL,NULL);
  cell_grp_id = H5Gopen(file_id,"Cells");
  H5Eset_auto(old_func,old_client_data);
  if (cell_grp_id < 0)      */ // does not matter since file is being overwritten
                                     // < 0 sets to default hint
    cell_grp_id = H5Gcreate(file_id,"Cells",-1) ;

// Write cells
  // Create data space 
  int rank = 2;
  hsize_t dims[3];
  hsize_t max_dims[3];
  dims[0] = grid->getNumberOfCells();
  dims[1] = 8;
  max_dims[0] = grid->getNumberOfCells();
  max_dims[1] = 8;
  data_space_id = H5Screate_simple(rank,dims,max_dims);

  // Create data set
  prop_id = H5Pcreate(H5P_DATASET_CREATE);
  // here is where you set chunking, shuffle, compression, deflate, fill_value, etc.
  data_set_id = H5Dcreate(cell_grp_id,"VertexIds",H5T_NATIVE_INT,data_space_id,
                          prop_id);
  H5Pclose(prop_id);

  prop_id = H5Pcreate(H5P_DATASET_XFER);
#ifdef PARALLEL
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
#endif

  int *vertex_ids = new int[grid->getNumberOfCells()*8];
  for (int i=0; i<grid->getNumberOfCells()*8; i++)
    vertex_ids[i] = -999;
  for (int icell=0; icell<grid->getNumberOfCells(); icell++)
    for (int ivert=0; ivert<grid->cells[icell].vertices[0]; ivert++)
      vertex_ids[ivert+icell*8] = grid->cells[icell].vertices[ivert+1]; // vertices[0] store the number of verts

  H5Dwrite(data_set_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,prop_id,vertex_ids);
  H5Pclose(prop_id);

  delete [] vertex_ids;

  H5Dclose(data_set_id);
  H5Sclose(data_space_id);
  H5Gclose(cell_grp_id);

// Open vertex group  
/*
  herr_t (*old_func)(void*);
  void  *old_client_data;
  H5Eget_auto(&old_func, &old_client_data);
  H5Eset_auto(NULL,NULL);
  cell_grp_id = H5Gopen(file_id,"Cells");
  H5Eset_auto(old_func,old_client_data);
  if (cell_grp_id < 0)      */ // does not matter since file is being overwritten
                                     // < 0 sets to default hint
    vert_grp_id = H5Gcreate(file_id,"Vertices",-1) ;

// Write vertices
  // Create data space 
  rank = 1;
  dims[0] = grid->getNumberOfVertices();
  max_dims[0] = grid->getNumberOfVertices();
  data_space_id = H5Screate_simple(rank,dims,max_dims);

// vertex natural ids
  // Create data set
  prop_id = H5Pcreate(H5P_DATASET_CREATE);
  // here is where you set chunking, shuffle, compression, deflate, fill_value, etc.
  data_set_id = H5Dcreate(vert_grp_id,"Ids",H5T_NATIVE_INT,data_space_id,
                          prop_id);
  H5Pclose(prop_id);

  prop_id = H5Pcreate(H5P_DATASET_XFER);
#ifdef PARALLEL
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
#endif

  vertex_ids = new int[grid->getNumberOfVertices()];
  for (int ivert=0; ivert<grid->getNumberOfVertices(); ivert++)
    vertex_ids[ivert] = -999;
  for (int ivert=0; ivert<grid->getNumberOfVertices(); ivert++)
    vertex_ids[ivert] = grid->vertices[ivert].getIdNatural();

  H5Dwrite(data_set_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,prop_id,vertex_ids);
  H5Pclose(prop_id);
  H5Dclose(data_set_id);

  delete [] vertex_ids;

// x-coordinate
  // Create data set
  prop_id = H5Pcreate(H5P_DATASET_CREATE);
  // here is where you set chunking, shuffle, compression, deflate, fill_value, etc.
  data_set_id = H5Dcreate(vert_grp_id,"X-Coordinates",H5T_NATIVE_DOUBLE,data_space_id,
                          prop_id);
  H5Pclose(prop_id);

  prop_id = H5Pcreate(H5P_DATASET_XFER);
#ifdef PARALLEL
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
#endif

  double *coordinate = new double[grid->getNumberOfVertices()];
  for (int ivert=0; ivert<grid->getNumberOfVertices(); ivert++)
    coordinate[ivert] = -999.;
  for (int ivert=0; ivert<grid->getNumberOfVertices(); ivert++)
    coordinate[ivert] = grid->vertices[ivert].getX();

  H5Dwrite(data_set_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,prop_id,coordinate);
  H5Pclose(prop_id);
  H5Dclose(data_set_id);

// y-coordinate
  // Create data set
  prop_id = H5Pcreate(H5P_DATASET_CREATE);
  // here is where you set chunking, shuffle, compression, deflate, fill_value, etc.
  data_set_id = H5Dcreate(vert_grp_id,"Y-Coordinates",H5T_NATIVE_DOUBLE,data_space_id,
                          prop_id);
  H5Pclose(prop_id);

  prop_id = H5Pcreate(H5P_DATASET_XFER);
#ifdef PARALLEL
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
#endif

  coordinate = new double[grid->getNumberOfVertices()];
  for (int ivert=0; ivert<grid->getNumberOfVertices(); ivert++)
    coordinate[ivert] = -999.;
  for (int ivert=0; ivert<grid->getNumberOfVertices(); ivert++)
    coordinate[ivert] = grid->vertices[ivert].getY();

  H5Dwrite(data_set_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,prop_id,coordinate);
  H5Pclose(prop_id);
  H5Dclose(data_set_id);

// z-coordinate
  // Create data set
  prop_id = H5Pcreate(H5P_DATASET_CREATE);
  // here is where you set chunking, shuffle, compression, deflate, fill_value, etc.
  data_set_id = H5Dcreate(vert_grp_id,"Z-Coordinates",H5T_NATIVE_DOUBLE,data_space_id,
                          prop_id);
  H5Pclose(prop_id);

  prop_id = H5Pcreate(H5P_DATASET_XFER);
#ifdef PARALLEL
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
#endif

  for (int ivert=0; ivert<grid->getNumberOfVertices(); ivert++)
    coordinate[ivert] = -999.;
  for (int ivert=0; ivert<grid->getNumberOfVertices(); ivert++)
    coordinate[ivert] = grid->vertices[ivert].getZ();

  H5Dwrite(data_set_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,prop_id,coordinate);
  H5Pclose(prop_id);
  H5Dclose(data_set_id);

  H5Sclose(data_space_id);
  H5Gclose(vert_grp_id);

  H5Fclose(file_id);
#endif
}

void Output::printHDFMesh2() {

  PetscPrintf(PETSC_COMM_WORLD,"Printing HDF grid.\n");

  HDF *file = new HDF("grid.h5",1);
  int compress = 1;

// write cells
  file->createGroup("Cells");

  file->createDataSpace(2,grid->getGlobalNumberOfCells(),8,0);
  file->createDataSet("VertexIds",H5T_NATIVE_INT,compress);

  int *vertex_ids = grid->getLocalCellVertexNaturalIDs();

  file->writeInt(vertex_ids);
  delete [] vertex_ids;

  file->closeDataSet();
  file->closeDataSpace();

// natural ids
  file->createDataSpace(1,grid->getGlobalNumberOfCells(),0,0);
  file->createDataSet("NaturalIds",H5T_NATIVE_INT,compress);

  int *natural_ids = grid->getLocalCellNaturalIDs();
  file->writeInt(natural_ids);
  delete [] natural_ids;

  file->closeDataSet();

// materials
  // use same data space
  file->createDataSet("Materials",H5T_NATIVE_INT,compress);

  int *material_ids = grid->getLocalCellMaterialNaturalIDs();

  file->writeInt(material_ids);
  delete [] material_ids;

  file->closeDataSet();
  file->closeDataSpace();

  file->closeGroup();

// vertices
  file->createGroup("Vertices");
  file->createDataSpace(1,grid->getGlobalNumberOfVertices(),0,0);

// natural ids
  file->createDataSet("NaturalIds",H5T_NATIVE_INT,compress);

  natural_ids = grid->getLocalCellVertexNaturalIDs();

  file->writeInt(natural_ids);
  file->closeDataSet();

  delete [] natural_ids;

// x-coordinate
  file->createDataSet("X-Coordinates",H5T_NATIVE_DOUBLE,compress);

  double *coordinate = new double[grid->getNumberOfVertices()];
  for (int ivert=0; ivert<grid->getNumberOfVertices(); ivert++)
    coordinate[ivert] = -999.;
  for (int ivert=0; ivert<grid->getNumberOfVertices(); ivert++)
    coordinate[ivert] = grid->vertices[ivert].getX();

  file->writeDouble(coordinate);
  file->closeDataSet();

// y-coordinate
  file->createDataSet("Y-Coordinates",H5T_NATIVE_DOUBLE,compress);

  coordinate = new double[grid->getNumberOfVertices()];
  for (int ivert=0; ivert<grid->getNumberOfVertices(); ivert++)
    coordinate[ivert] = -999.;
  for (int ivert=0; ivert<grid->getNumberOfVertices(); ivert++)
    coordinate[ivert] = grid->vertices[ivert].getY();

  file->writeDouble(coordinate);
  file->closeDataSet();

// z-coordinate
  file->createDataSet("Z-Coordinates",H5T_NATIVE_DOUBLE,compress);

  coordinate = new double[grid->getNumberOfVertices()];
  for (int ivert=0; ivert<grid->getNumberOfVertices(); ivert++)
    coordinate[ivert] = -999.;
  for (int ivert=0; ivert<grid->getNumberOfVertices(); ivert++)
    coordinate[ivert] = grid->vertices[ivert].getZ();

  file->writeDouble(coordinate);
  file->closeDataSet();

  delete [] coordinate;

  file->closeDataSpace();
  file->closeGroup();

// boundary connections
  file->createGroup("BoundaryConnections");

  BoundarySet *cur_set = grid->boundary_sets;
  while (cur_set) {
    file->createGroup(cur_set->name);
// condition id
    char *name = cur_set->condition->getName();
    file->writeString("Condition",name);
// face vertex ids
    file->createDataSpace(2,cur_set->getNumberOfConnections(),4,0);
    file->createDataSet("FaceVertexIds",H5T_NATIVE_INT,compress);

    int *vertex_ids = new int[4*cur_set->getNumberOfConnections()];
    Connection *cur_conn = cur_set->list;
    int count = 0;
    int face_vertices[4];

    while (cur_conn) {
      grid->cells[cur_conn->cell].getHexFaceVertices(cur_conn->face,face_vertices);
      for (int i=0; i<4; i++)
        vertex_ids[count*4+i] = face_vertices[i];
      count++;
      cur_conn = cur_conn->next;
    }
    file->writeInt(vertex_ids);

    file->closeDataSet();
    file->closeDataSpace();
    file->closeGroup();
    cur_set = cur_set->next;
  }
  file->closeGroup();

// conditions
  file->createGroup("Conditions");

  Condition *cur_conn = Condition::list;
  while (cur_conn) {

    char *name = cur_conn->getName();
    file->createGroup(name);

// times
    file->createDataSpace(1,cur_conn->max_time_index,0,0);
    file->createDataSet("Times",H5T_NATIVE_DOUBLE,compress);

    file->writeDouble(cur_conn->times);

    file->closeDataSet();

// scalar values
// same data space
    file->createDataSet("Values",H5T_NATIVE_DOUBLE,compress);

    file->writeDouble(cur_conn->scalars);

    file->closeDataSet();

    file->closeDataSpace();
    file->closeGroup();
    cur_conn = cur_conn->next;
  } 

  file->closeGroup();

  delete file;

}

Output::~Output() {

}
