//
// $Id: test_gradient.C 1832 2005-08-23 21:10:00Z bphilip $
// $Revision: 1832 $
// $Date: 2005-08-23 15:10:00 -0600 (Tue, 23 Aug 2005) $
//

#include "SAMRAI_config.h"

#include <iostream>

#include <fstream>
using namespace std;

#include <sys/stat.h>

/*
 * SAMRAI headers.
 */
#include "BoundaryBox.h"
#include "BoxArray.h"
#include "BoxList.h"
#include "BergerRigoutsos.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellVariable.h"
#include "FaceData.h"
#include "FaceVariable.h"
#include "GriddingAlgorithm.h"
#include "Index.h"
#include "tbox/InputDatabase.h"
#include "tbox/InputManager.h"
#include "IntVector.h"
#include "LoadBalancer.h"
#include "tbox/MPI.h"
#include "Patch.h"
#include "PatchGeometry.h"
#include "tbox/PIO.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "tbox/SAMRAIManager.h"
#include "StandardTagAndInitialize.h"
#include "tbox/Utilities.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "PETSc_SAMRAIVectorReal.h"
extern "C"{
#include "assert.h"
}
/*
 * Application header.
 */
#include "BoundaryConditionStrategy.h"
#include "BogusTagAndInitStrategy.h"
#include "fc_interface.h"
#include "SAMRAIDriver.h"
/*#include "pims_local_struct.h"*/
/*
 * Ghost cell width for variables that need them.
 */

#define GHOST_CELL_WIDTH (1)

int main( int argc, char *argv[] ) 
{
   double max_error = 0.0;
   double l2_error = 0.0;
   string input_file;
   string log_file;
   
   tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy;

   tbox::MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();
   int ierr = PetscInitializeNoArguments();
   PetscInitializeFortran();

   /*
    * Process command line arguments and dump to log file.
    */
   processCommandLine(argc, argv, input_file, log_file);

   tbox::PIO::logOnlyNodeZero(log_file);

   /*
    * Create input database and parse all data in input file.  This
    * parsing allows us to subsequently extract individual sections.
    */
   tbox::Pointer<tbox::Database> input_db = new tbox::InputDatabase("input_db");
   tbox::InputManager::getManager()->parseInputFile(input_file, input_db);

   /*
    * Create an application object.  Some TagAndInitStrategy must be
    * provided in order to build an object that specifies cells that
    * need refinement.  Here an empty object is provided, since a
    * prescribed set of refinement regions are read in from the input
    * file; it would be a useful exercise to fill in the
    * applyGradientDetector method and to generate custom refinement
    * regions.
    */
   BogusTagAndInitStrategy* test_object = new BogusTagAndInitStrategy();

   /*
    * Create the AMR hierarchy and initialize it
    */
   initializeAMRHierarchy(input_db,
			  test_object,
			  hierarchy);
 

   void *p_pflowhierarchy=NULL;
   f_create_hierarchy_data_(&p_pflowhierarchy);

   void *p_integrator = NULL;
   f_create_integrator_(&p_integrator);

   void *p_samr_hierarchy=(void *)(hierarchy.getPointer());
   gridparameters *params = new gridparameters;
   params->p_grid = p_pflowhierarchy;
   params->p_timestep = p_integrator;
   params->igeom = 1;
   params->nx = 1;
   params->ny = 1;
   params->nz = 1;
   params->npx = 1;
   params->npy = 1;
   params->npz = 1;
   params->nphase = 1;
   params->nlevels = hierarchy->getNumberLevels();
   params->usesamrai = PETSC_TRUE;
   params->p_samr_hierarchy = p_samr_hierarchy;

   f_initialize_hierarchy_data_((void **)&params);

   f_setup_hierarchy_data_((void **)&p_pflowhierarchy, &p_integrator);

   /* 
    * That's all, folks!
    */
   PetscFinalize();
   tbox::SAMRAIManager::shutdown();
   tbox::MPI::finalize();

   return(0);
}

/*
************************************************************************
*                                                                      *
*  Parse command line arguments, returning name of input file and log  *
*  file.                                                               *
*                                                                      *
************************************************************************
*/
void processCommandLine(int argc, 
                        char *argv[], 
                        string& input_file, 
                        string& log_file)
{
  if ( (argc != 3) ) {
    tbox::pout << "USAGE:  " << argv[0] << " <input file> <log file> " << endl;
    exit(-1);
  } else {
    input_file = argv[1];
    log_file = argv[2];
  }

  return;
}

/*
************************************************************************
*                                                                      *
* Generate a patch hierarchy from information specified in input       *
* database.                                                            *
*                                                                      *
************************************************************************
*/
void initializeAMRHierarchy(tbox::Pointer<tbox::Database> &input_db,
			    mesh::StandardTagAndInitStrategy<NDIM>* user_tagging_strategy,
			    tbox::Pointer<hier::PatchHierarchy<NDIM> > &hierarchy)
{
   /*
    * Create geometry object.  This specifies the index space of the
    * coarsest level, as well as its physical (Cartesian) coordinates.
    */
   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry = 
      new geom::CartesianGridGeometry<NDIM>("CartesianGeometry",
                                input_db->getDatabase("CartesianGeometry"));

   /*
    * Create patch hierarchy.
    */
   hierarchy = new hier::PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);

   /* 
    * A mesh::GriddingAlgorithm<NDIM> is used to build the initial grid hierarchy.
    * Classes for tagging cells that need refinement, generation of
    * boxes from these tagged cells, and load balancing the grid
    * hierarchy are needed to build the mesh::GriddingAlgorithm<NDIM>.
    *
    * First build the object used to tag cells that need refinement.
    */
    tbox::Pointer<mesh::StandardTagAndInitialize<NDIM> > error_detector = 
	new mesh::StandardTagAndInitialize<NDIM>( 
	    "CellTaggingMethod", 
	    user_tagging_strategy, 
	    input_db->getDatabase("StandardTagAndInitialize"));
    
   /*
    * Next, specify the built-in Berger-Rigoutsos method for
    * generating boxes from the tagged cells.
    */
   tbox::Pointer<mesh::BergerRigoutsos<NDIM> > box_generator = new mesh::BergerRigoutsos<NDIM>(); 

   /*
    * Next, specify the built-in uniform load balancer to distribute
    * patches across processors.
    */
   tbox::Pointer<mesh::LoadBalancer<NDIM> > load_balancer =
      new mesh::LoadBalancer<NDIM>(input_db->getDatabase("LoadBalancer"));

   /*
    * Finally, build the grid generator, registering the above
    * strategies for tagging cells, generating boxes, and load
    * balancing the calculation.
    */
   tbox::Pointer<mesh::GriddingAlgorithm<NDIM> > gridding_algorithm =
      new mesh::GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                            input_db->getDatabase("GriddingAlgorithm"),
                            error_detector,
                            box_generator,
                            load_balancer);

   /*
    * Build an initial grid hierarchy.  Note that in this simple
    * example we do not buffer the refinement regions.
    */
   gridding_algorithm->makeCoarsestLevel(hierarchy, 0.0);
   
   bool done = false;
   bool initial_time = true;
   for (int ln = 0;
        gridding_algorithm->levelCanBeRefined(ln) && !done; 
        ln++) {
       gridding_algorithm->makeFinerLevel(hierarchy,
                                          0.0,
                                          initial_time,
                                          0);
       done = !(hierarchy->finerLevelExists(ln));
   }
}

