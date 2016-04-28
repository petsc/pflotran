#include "Output.h"

Output::Output(Grid *g) {

  grid = g;

}

void Output::printGMSGrid() {

  FILE *fp = NULL;
  
  char filename[128];
  char string[128];
  
  grid->getFilenamePrefix(string);
  strcpy(filename,"gms_");
  strcat(filename,string);
  strcat(filename,".3dg");

  PetscPrintf(PETSC_COMM_WORLD,"Printing GMS grid.\n");
  PetscFOpen(PETSC_COMM_WORLD,filename,"w",&fp);
  PetscFPrintf(PETSC_COMM_WORLD,fp,"GRID3D\n");
  PetscFPrintf(PETSC_COMM_WORLD,fp,"TYPE 1\n");
  PetscFPrintf(PETSC_COMM_WORLD,fp,"IJK +y +x +z\n");
  PetscFPrintf(PETSC_COMM_WORLD,fp,"NUMBERING 0\n");
  PetscReal *ptr = grid->getOriginPtr();
  PetscFPrintf(PETSC_COMM_WORLD,fp,"ORIGIN %.8e %.8e %.8e\n",ptr[0],ptr[1],ptr[2]);
  ptr = NULL;
  PetscFPrintf(PETSC_COMM_WORLD,fp,"ROTZ %f\n",grid->getRotationDegrees());
  PetscFPrintf(PETSC_COMM_WORLD,fp,"DIM %d %d %d\n",grid->getNx()+1,
               grid->getNy()+1,grid->getNz()+1);

  PetscReal sum = 0.;
  PetscFPrintf(PETSC_COMM_WORLD,fp,"%.8e\n",sum);
  for (PetscInt i=0; i<grid->getNx(); i++) {
    sum += grid->getDx(i);
    PetscFPrintf(PETSC_COMM_WORLD,fp,"%.8e\n",sum);
  }
  sum = 0.;
  PetscFPrintf(PETSC_COMM_WORLD,fp,"%.8e\n",sum);
  for (PetscInt j=0; j<grid->getNy(); j++) {
    sum += grid->getDy(j);
    PetscFPrintf(PETSC_COMM_WORLD,fp,"%.8e\n",sum);
  }
  sum = 0.;
  PetscFPrintf(PETSC_COMM_WORLD,fp,"%.8e\n",sum);
  for (PetscInt k=0; k<grid->getNz(); k++) {
    sum += grid->getDz(k);
    PetscFPrintf(PETSC_COMM_WORLD,fp,"%.8e\n",sum);
  }
  PetscFPrintf(PETSC_COMM_WORLD,fp,"MAT\n");
  Vec v = grid->getGridCellMaterialIDs();
  writeIntVectorInNaturalOrder(fp,v,0);
  PetscFPrintf(PETSC_COMM_WORLD,fp,"ACTIVE\n");
  v = grid->getGridCellActivities();
  writeIntVectorInNaturalOrder(fp,v,0);
  VecDestroy(&v);

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
  PetscInt i;
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

void Output::writeIntVectorInNaturalOrder(FILE *fp, Vec v, PetscInt one_per_line) {
  Vec natural;
  PetscInt num_per_line = 20;
  PetscReal *v_ptr = NULL;
  grid->getVectorNatural(&natural);
  grid->globalToNatural(v,natural);
  PetscInt max_cells;
  MPI_Allreduce(&grid->num_cells_local,&max_cells,1,MPIU_INT,MPI_MAX,PETSC_COMM_WORLD);
  PetscInt *values = new PetscInt[max_cells];
  VecGetArray(natural,&v_ptr);
  for (PetscInt i=0; i<grid->num_cells_local; i++) 
    values[i] = PetscInt(v_ptr[i]);
  VecRestoreArray(natural,&v_ptr);
  VecDestroy(&natural);
  if (myrank == 0) {
    PetscInt count = 0;
    for (PetscInt i=0; i<grid->num_cells_local; i++) {
      if (one_per_line) {
        PetscFPrintf(PETSC_COMM_WORLD,fp,"%d\n",values[i]);
      }
      else {
        PetscFPrintf(PETSC_COMM_WORLD,fp,"%d ",values[i]);
        if (++count%num_per_line == 0) PetscFPrintf(PETSC_COMM_WORLD,fp,"\n");
      }
    }
    for (PetscInt iproc=1; iproc<commsize; iproc++) {
      MPI_Status status;
      MPI_Probe(iproc,MPI_ANY_TAG,PETSC_COMM_WORLD,&status);
      PetscInt nrecv = status.MPI_TAG;
      MPI_Recv(values,nrecv,MPIU_INT,iproc,MPI_ANY_TAG,PETSC_COMM_WORLD,&status);
      for (PetscInt i=0; i<nrecv; i++) {
        if (one_per_line) {
          PetscFPrintf(PETSC_COMM_WORLD,fp,"%d\n",values[i]);
        }
        else {
          PetscFPrintf(PETSC_COMM_WORLD,fp,"%d ",values[i]);
          if (++count%num_per_line == 0) PetscFPrintf(PETSC_COMM_WORLD,fp,"\n");
        }
      }
    }
    if (!one_per_line && count%num_per_line != 0) PetscFPrintf(PETSC_COMM_WORLD,fp,"\n");
  }
  else {
    MPI_Send(values,grid->num_cells_local,MPIU_INT,0,grid->num_cells_local,
             PETSC_COMM_WORLD);
  }
  delete [] values;
}

void Output::printBoundarySets() {
  Vec v;
  PetscReal *v_ptr = NULL;
  grid->getVectorGlobal(&v);
  BoundarySet *cur_set = grid->boundary_sets;  
  while (cur_set) {
    VecGetArray(v,&v_ptr);
    for (PetscInt i=0; i<grid->num_cells_local; i++)
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
  VecDestroy(&v);
}

void Output::printHDFMaterialsAndRegions() {

#if defined(PETSC_HAVE_HDF5)
  PetscPrintf(PETSC_COMM_WORLD,"Printing HDF Materials and Regions.\n");

  PetscBool option_found;
  char filename[128];
  char string[128];
  
  grid->getFilenamePrefix(string);
  strcpy(filename,"input_");
  strcat(filename,string);
  strcat(filename,".h5");
  PetscOptionsGetString(PETSC_NULL,"-mdtout",filename,1024,&option_found);

  HDF *file = new HDF(filename,1);
  PetscInt compress = 0;

// materials
  PetscPrintf(PETSC_COMM_WORLD,"Materials\n");
  file->createGroup("Materials");

  file->createFileSpace(1,grid->getNumberOfCellsGlobal(),NULL,NULL);
  PetscPrintf(PETSC_COMM_WORLD,"Cell Ids\n");
  file->createDataSet("Cell Ids",HDF_NATIVE_INT,compress);

  PetscInt *cell_ids = grid->getCellIdsNatural1Based();
  grid->convertLocalCellDataGtoN(cell_ids);

  file->setHyperSlab(grid->getNumberOfCellsLocal());
  file->createMemorySpace(1,grid->getNumberOfCellsLocal(),NULL,NULL);
  file->writeInt(cell_ids);

  delete [] cell_ids;
  cell_ids = NULL;

  file->closeDataSet();

  // use same data space
  PetscPrintf(PETSC_COMM_WORLD,"Material Ids\n");
  file->createDataSet("Material Ids",HDF_NATIVE_INT,compress);

  PetscInt *material_ids = grid->getCellMaterialIds();
  grid->convertLocalCellDataGtoN(material_ids);

  file->setHyperSlab(grid->getNumberOfCellsLocal());
  file->createMemorySpace(1,grid->getNumberOfCellsLocal(),NULL,NULL);
  file->writeInt(material_ids);
  delete [] material_ids;
  material_ids = NULL;

  file->closeDataSet();
  file->closeDataSpaces();
  file->closeGroup();

// boundary connections
  PetscPrintf(PETSC_COMM_WORLD,"Regions\n");
  file->createGroup("Regions");

  BoundarySet *cur_set = grid->boundary_sets;
  while (cur_set) {
    PetscPrintf(PETSC_COMM_WORLD," %s\n",cur_set->name);
//    cur_set->printInfo();
// cell ids
    PetscInt num_connections_global = cur_set->getNumberOfConnectionsGlobal();
    PetscInt num_connections_local = cur_set->getNumberOfConnectionsLocal();

    if (num_connections_global == 0) {
      PetscPrintf(PETSC_COMM_WORLD,"   Boundary Condition: %s has not connections, skipping.\n",cur_set->name);
    }
    else {
      
      file->createGroup(cur_set->name);

  // hdf5 cannot handle a zero-length memory space and hyperslab.  Therefore,
  // trick it by setting then length to 1 and writing independently (non-
  // collective) for those procs with num_connections_local > 0).
      PetscInt trick_hdf5_local = num_connections_local > 0 ? 
                                               num_connections_local : 1;
      PetscInt trick_hdf5_global = num_connections_global > 0 ? 
                                               num_connections_global : 1;

      file->createDataSpace(1,trick_hdf5_global,0,0);
      PetscPrintf(PETSC_COMM_WORLD,"  Cell Ids\n");
      file->createDataSet("Cell Ids",HDF_NATIVE_INT,compress);

      if (num_connections_global > 0) {

        PetscInt *cell_ids = cur_set->getCellIdsLocal1Based();
      // convert cell_ids in local numbering to natural, no need to reorder for now
      // since 1-based, need to subtract 1 then add 1
        for (PetscInt i=0; i<num_connections_local; i++) {
          cell_ids[i] = grid->cell_mapping_ghosted_to_natural[
                              grid->cell_mapping_local_to_ghosted[cell_ids[i]-1]]+1;
        }

        file->setHyperSlab(num_connections_local);
        file->createMemorySpace(1,num_connections_local,NULL,NULL);
        if (num_connections_local > 0) {
          file->writeInt(cell_ids,INDEPENDENT);
        }
        delete [] cell_ids;
        cell_ids = NULL;
      }
      file->closeDataSet();
      file->closeDataSpaces();

      file->createDataSpace(1,trick_hdf5_global,0,0);
      PetscPrintf(PETSC_COMM_WORLD,"  Face Ids\n");
      file->createDataSet("Face Ids",HDF_NATIVE_INT,compress);

      if (num_connections_global > 0) {

        PetscInt *face_ids = cur_set->getFaceIdsLocal();

        file->setHyperSlab(num_connections_local);
        file->createMemorySpace(1,num_connections_local,NULL,NULL);
        if (num_connections_local > 0) {
          file->writeInt(face_ids,INDEPENDENT);
        }
        delete [] face_ids;
        face_ids = NULL;
      }
      file->closeDataSet();
      file->closeDataSpaces();

      file->closeGroup();
    }
    cur_set = cur_set->next;
  }
  file->closeGroup();

  delete file;
  PetscPrintf(PETSC_COMM_WORLD,"Done with HDF5 Materials and Regions\n");
#endif

}

void Output::printHDFMesh() {

#if defined(PETSC_HAVE_HDF5)
  PetscPrintf(PETSC_COMM_WORLD,"Printing HDF grid.\n");

  HDF *file = new HDF("grid.h5",1);
  PetscInt compress = 0;

// write cells
  PetscPrintf(PETSC_COMM_WORLD,"Cells\n");
  file->createGroup("Cells");

  // just to be sure everything is closed.
  file->closeDataSpaces();

  file->createFileSpace(1,grid->getNumberOfCellsGlobal(),NULL,NULL);
  PetscPrintf(PETSC_COMM_WORLD,"CellIds\n");
  file->createDataSet("CellIds",HDF_NATIVE_INT,compress);

  PetscInt *cell_ids = grid->getCellIdsNatural();
  grid->convertLocalCellDataGtoN(cell_ids);

  file->setHyperSlab(grid->getNumberOfCellsLocal());
  file->createMemorySpace(1,grid->getNumberOfCellsLocal(),NULL,NULL);
  file->writeInt(cell_ids);

  delete [] cell_ids;
  cell_ids = NULL;

  file->closeDataSet();

// natural ids
  file->createFileSpace(1,grid->getNumberOfCellsGlobal(),NULL,NULL);
  PetscPrintf(PETSC_COMM_WORLD,"NaturalIds\n");
  file->createDataSet("NaturalIds",HDF_NATIVE_INT,compress);

  PetscInt *natural_ids = grid->getCellIdsNatural();
  grid->convertLocalCellDataGtoN(natural_ids);

  file->setHyperSlab(grid->getNumberOfCellsLocal());
  file->createMemorySpace(1,grid->getNumberOfCellsLocal(),NULL,NULL);
  file->writeInt(natural_ids);

  delete [] natural_ids;
  natural_ids = NULL;

  file->closeDataSet();

// materials
  // use same data space
  PetscPrintf(PETSC_COMM_WORLD,"Materials\n");
  file->createDataSet("Materials",HDF_NATIVE_INT,compress);

  PetscInt *material_ids = grid->getCellMaterialIds();
  grid->convertLocalCellDataGtoN(material_ids);

  file->setHyperSlab(grid->getNumberOfCellsLocal());
  file->createMemorySpace(1,grid->getNumberOfCellsLocal(),NULL,NULL);
  file->writeInt(material_ids);
  delete [] material_ids;
  material_ids = NULL;

  file->closeDataSet();
  file->closeDataSpaces();

  // cell vertices
  PetscInt num_cells_global = grid->getNumberOfCellsGlobal();
  file->createFileSpace(2,num_cells_global,8,NULL);
  PetscPrintf(PETSC_COMM_WORLD,"CellVertices\n");
  file->createDataSet("CellVertices",HDF_NATIVE_INT,compress);
  file->createMemorySpace(1,grid->getNumberOfCellsGlobal(),NULL,NULL);

  PetscInt offset = 0;
  PetscInt num_cells_local = grid->getNumberOfCellsLocal();
  MPI_Exscan(&num_cells_local,&offset,1,MPIU_INT,MPI_SUM,PETSC_COMM_WORLD);
  if (offset > num_cells_global)
    printf("Proc[%d]: ERROR - offset(%d) > num_cells_global(%d)\n",
           myrank, offset, num_cells_global);
  if (offset == num_cells_global) {
    printf("Proc[%d]: WARNING - offset(%d) == num_cells_global(%d)\n",
           myrank, offset, num_cells_global);
    offset--;
  }

  for (PetscInt ivert=0; ivert<8; ivert++) {

    PetscInt *vertex_ids = grid->getCellVertexIds(ivert+1);
    grid->convertLocalCellDataGtoN(vertex_ids);

    PetscInt start[3] = {offset,ivert,0};
    PetscInt stride[3] = {1,8,1};
    PetscInt count[3] = {num_cells_local,1,1};
    PetscInt block[3] = {1,1,1};
    file->setHyperSlab(start,stride,count,NULL);
    file->createMemorySpace(1,count[0]*count[1]*count[2],NULL,NULL);
    file->writeInt(vertex_ids);
    delete [] vertex_ids;
    vertex_ids = NULL;
  }

  file->closeDataSet();
  file->closeDataSpaces();

  // close cell group
  file->closeGroup();

// vertices
  PetscPrintf(PETSC_COMM_WORLD,"Vertices\n");
  file->createGroup("Vertices");

  file->createDataSpace(1,grid->getNumberOfVerticesGlobal(),0,0);

// natural ids

  PetscPrintf(PETSC_COMM_WORLD,"NaturalIds\n");
  file->createDataSet("NaturalIds",HDF_NATIVE_INT,compress); 

  PetscInt num_print_vertices_local = grid->getVertexIdsNaturalLocal(&natural_ids);

  file->setHyperSlab(num_print_vertices_local);
  file->createMemorySpace(1,num_print_vertices_local,NULL,NULL);
  file->writeInt(natural_ids);

  delete [] natural_ids;
  natural_ids = NULL;

  file->closeDataSet();

// x-coordinate
  PetscPrintf(PETSC_COMM_WORLD,"X-Coordinates\n");
  file->createDataSet("X-Coordinates",H5T_NATIVE_DOUBLE,compress);

  PetscReal *coordinates = NULL;
  num_print_vertices_local = 
                grid->getVertexCoordinatesNaturalLocal(&coordinates,0); // 0 = X

  file->setHyperSlab(num_print_vertices_local);
  file->createMemorySpace(1,num_print_vertices_local,NULL,NULL);
  file->writeDouble(coordinates);
  file->closeDataSet();

  delete [] coordinates;
  coordinates = NULL;

// y-coordinate
  PetscPrintf(PETSC_COMM_WORLD,"Y-Coordinates\n");
  file->createDataSet("Y-Coordinates",H5T_NATIVE_DOUBLE,compress);

  num_print_vertices_local = 
                grid->getVertexCoordinatesNaturalLocal(&coordinates,1); // 1 = Y

  file->setHyperSlab(num_print_vertices_local);
  file->createMemorySpace(1,num_print_vertices_local,NULL,NULL);
  file->writeDouble(coordinates);
  file->closeDataSet();

  delete [] coordinates;
  coordinates = NULL;

// z-coordinate
  PetscPrintf(PETSC_COMM_WORLD,"Z-Coordinates\n");
  file->createDataSet("Z-Coordinates",H5T_NATIVE_DOUBLE,compress);

  num_print_vertices_local = 
                grid->getVertexCoordinatesNaturalLocal(&coordinates,2); // 2 = Z

  file->setHyperSlab(num_print_vertices_local);
  file->createMemorySpace(1,num_print_vertices_local,NULL,NULL);
  file->writeDouble(coordinates);
  file->closeDataSet();

  delete [] coordinates;
  coordinates = NULL;

  file->closeDataSpaces();

  file->closeGroup();

// boundary connections
  PetscPrintf(PETSC_COMM_WORLD,"BoundaryConnections\n");
  file->createGroup("BoundaryConnections");

  BoundarySet *cur_set = grid->boundary_sets;
  while (cur_set) {
    PetscPrintf(PETSC_COMM_WORLD," %s\n",cur_set->name);
    file->createGroup(cur_set->name);
// condition id
    if (cur_set->condition) {
      char *name = cur_set->condition->getName();
      PetscPrintf(PETSC_COMM_WORLD,"  Condition %s\n",name);
      file->writeString("Condition",name,INDEPENDENT);
    }

// cell ids
    PetscInt num_connections_global = cur_set->getNumberOfConnectionsGlobal();
    PetscInt num_connections_local = cur_set->getNumberOfConnectionsLocal();
// hdf5 cannot handle a zero-length memory space and hyperslab.  Therefore,
// trick it by setting then length to 1 and writing independently (non-
// collective) for those procs with num_connections_local > 0).
    PetscInt trick_hdf5_local = num_connections_local > 0 ? 
                                             num_connections_local : 1;
    PetscInt trick_hdf5_global = num_connections_global > 0 ? 
                                             num_connections_global : 1;

    file->createDataSpace(1,trick_hdf5_global,0,0);
    PetscPrintf(PETSC_COMM_WORLD,"  CellIds\n");
    file->createDataSet("CellIds",HDF_NATIVE_INT,compress);

    if (num_connections_global > 0) {

      cell_ids = cur_set->getCellIdsLocal();
    // convert cell_ids in local numbering to natural, no need to reorder for now
      for (PetscInt i=0; i<num_connections_local; i++) {
        cell_ids[i] = grid->cell_mapping_ghosted_to_natural[
                            grid->cell_mapping_local_to_ghosted[cell_ids[i]]];
      }

      file->setHyperSlab(num_connections_local);
      file->createMemorySpace(1,num_connections_local,NULL,NULL);
      if (num_connections_local > 0) {
        file->writeInt(cell_ids,INDEPENDENT);
      }
      delete [] cell_ids;
      cell_ids = NULL;
    }
    file->closeDataSet();
    file->closeDataSpaces();

    // face vertex ids

    file->createFileSpace(2,trick_hdf5_global,4,NULL);
    PetscPrintf(PETSC_COMM_WORLD,"  FaceVertexIds\n");
    file->createDataSet("FaceVertexIds",HDF_NATIVE_INT,compress);
    if (num_connections_global > 0) {
      PetscInt offset = 0;
      MPI_Exscan(&num_connections_local,&offset,1,MPIU_INT,MPI_SUM,
                 PETSC_COMM_WORLD);
      if (offset > num_connections_global)
        printf("Proc[%d]: ERROR - offset(%d) > num_connections_global(%d)\n",
               myrank, offset, num_connections_global);
      if (offset == num_connections_global) {
        printf("Proc[%d]: WARNING - offset(%d) == num_connections_global(%d)\n",
               myrank, offset, num_connections_global);
        offset--;
      }
      for (PetscInt ivert=0; ivert<4; ivert++) {

        PetscInt *vertex_ids = cur_set->getFaceVertexIds(ivert);
        PetscInt start[3] = {offset,ivert,0};
        PetscInt stride[3] = {1,4,1};
        PetscInt count[3] = {trick_hdf5_local,1,1};
        PetscInt block[3] = {1,1,1};
        file->setHyperSlab(start,stride,count,NULL);
        file->createMemorySpace(1,count[0]*count[1]*count[2],NULL,NULL);
        if (num_connections_local > 0) {
          file->writeInt(vertex_ids,INDEPENDENT);
        }
        delete [] vertex_ids;
        vertex_ids = NULL;
      }
    }

    file->closeDataSet();
    file->closeDataSpaces();
    file->closeGroup();
    cur_set = cur_set->next;
  }
  file->closeGroup();

// conditions
  PetscPrintf(PETSC_COMM_WORLD,"Conditions\n");
  file->createGroup("Conditions");

  Condition *cur_conn = Condition::list;
  while (cur_conn) {

    PetscPrintf(PETSC_COMM_WORLD," %s\n",cur_conn->getName());
    char *name = cur_conn->getName();
    file->createGroup(name);

// times
    file->createDataSpace(1,cur_conn->max_time_index,0,0);
    PetscPrintf(PETSC_COMM_WORLD," Times\n");
    file->createDataSet("Times",H5T_NATIVE_DOUBLE,compress);
 
    if (myrank == 0) file->writeDouble(cur_conn->times,INDEPENDENT);

    file->closeDataSet();

// scalar values
// same data space
    PetscPrintf(PETSC_COMM_WORLD," Values\n");
    file->createDataSet("Values",H5T_NATIVE_DOUBLE,compress);

    if (myrank == 0) file->writeDouble(cur_conn->scalars,INDEPENDENT);

    file->closeDataSet();

    file->closeDataSpaces();
    file->closeGroup();
    cur_conn = cur_conn->next;
  }

  file->closeGroup();

  delete file;
  PetscPrintf(PETSC_COMM_WORLD,"Done with HDF5 Mesh\n");
#endif

}

void Output::printHDFSieveMesh() {

#if defined(PETSC_HAVE_HDF5)
  PetscPrintf(PETSC_COMM_WORLD,"Printing HDF Sieve grid.\n");

  HDF *file = new HDF("sieve.h5",1);
  PetscInt compress = 0;

// write cells
  PetscPrintf(PETSC_COMM_WORLD,"Group: Connectivity\n");
  file->createGroup("Connectivity");

  // just to be sure everything is closed.
  file->closeDataSpaces();

  PetscInt num_cells_global = grid->getNumberOfCellsGlobal();

  file->writeAttribute("Number of Elements", num_cells_global);
  file->writeAttribute("Element Type", "hex");

  // cell vertices
  file->createFileSpace(2,num_cells_global,8,NULL);
  PetscPrintf(PETSC_COMM_WORLD,"Connectivity\n");
  file->createDataSet("Connectivity",HDF_NATIVE_INT,compress);
  file->createMemorySpace(1,grid->getNumberOfCellsGlobal(),NULL,NULL);

  PetscInt offset = 0;
  PetscInt num_cells_local = grid->getNumberOfCellsLocal();
  MPI_Exscan(&grid->num_cells_local,&offset,1,MPIU_INT,MPI_SUM,PETSC_COMM_WORLD);
  if (offset > num_cells_global)
    printf("Proc[%d]: ERROR - offset(%d) > num_cells_global(%d)\n",
           myrank, offset, num_cells_global);
  if (offset == num_cells_global) {
    printf("Proc[%d]: WARNING - offset(%d) == num_cells_global(%d)\n",
           myrank, offset, num_cells_global);
    offset--;
  }

  for (PetscInt ivert=0; ivert<8; ivert++) {

    PetscInt *vertex_ids = grid->getCellVertexIds(ivert+1);
    grid->convertLocalCellDataGtoN(vertex_ids);

    PetscInt start[3] = {offset,ivert,0};
    PetscInt stride[3] = {1,8,1};
    PetscInt count[3] = {num_cells_local,1,1};
    PetscInt block[3] = {1,1,1};
    file->setHyperSlab(start,stride,count,NULL);
    file->createMemorySpace(1,count[0]*count[1]*count[2],NULL,NULL);
    file->writeInt(vertex_ids);
    delete [] vertex_ids;
    vertex_ids = NULL;
  }

  file->closeDataSet();
  file->closeDataSpaces();

  // close cell group
  file->closeGroup();

// vertices
  PetscPrintf(PETSC_COMM_WORLD,"Group: Coordinates\n");
  file->createGroup("Coordinates");

  PetscInt num_vertices_global = grid->getNumberOfVerticesGlobal();
  file->writeAttribute("Number of Vertices", num_vertices_global);

  file->createDataSpace(1,num_vertices_global,0,0);

// x-coordinate
  PetscPrintf(PETSC_COMM_WORLD,"X-Coordinates\n");
  file->createDataSet("X-Coordinates",H5T_NATIVE_DOUBLE,compress);

  PetscReal *coordinates = NULL;
  PetscInt num_print_vertices_local = 
                grid->getVertexCoordinatesNaturalLocal(&coordinates,0); // 0 = X

  file->setHyperSlab(num_print_vertices_local);
  file->createMemorySpace(1,num_print_vertices_local,NULL,NULL);
  file->writeDouble(coordinates);
  file->closeDataSet();

  delete [] coordinates;
  coordinates = NULL;

// y-coordinate
  PetscPrintf(PETSC_COMM_WORLD,"Y-Coordinates\n");
  file->createDataSet("Y-Coordinates",H5T_NATIVE_DOUBLE,compress);

  num_print_vertices_local = 
                grid->getVertexCoordinatesNaturalLocal(&coordinates,1); // 1 = Y

  file->setHyperSlab(num_print_vertices_local);
  file->createMemorySpace(1,num_print_vertices_local,NULL,NULL);
  file->writeDouble(coordinates);
  file->closeDataSet();

  delete [] coordinates;
  coordinates = NULL;

// z-coordinate
  PetscPrintf(PETSC_COMM_WORLD,"Z-Coordinates\n");
  file->createDataSet("Z-Coordinates",H5T_NATIVE_DOUBLE,compress);

  num_print_vertices_local = 
                grid->getVertexCoordinatesNaturalLocal(&coordinates,2); // 2 = Z

  file->setHyperSlab(num_print_vertices_local);
  file->createMemorySpace(1,num_print_vertices_local,NULL,NULL);
  file->writeDouble(coordinates);
  file->closeDataSet();

  delete [] coordinates;
  coordinates = NULL;

  file->closeDataSpaces();

  file->closeGroup();

  delete file;
  PetscPrintf(PETSC_COMM_WORLD,"Done with HDF5 Sieve Mesh\n");
#endif

}

Output::~Output() {

}
