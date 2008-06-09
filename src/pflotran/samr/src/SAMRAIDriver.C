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
#include "tbox/SAMRAI_MPI.h"
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

#if 0
#include "PIMSApplication.h"
#include "PIMSApplicationParameters.h"
#include "TimeIntegratorParameters.h"
#include "ForwardEulerTimeIntegrator.h"
#endif

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

   tbox::SAMRAI_MPI::init(&argc, &argv);
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
 
   void *p_samr_hierarchy=(void *)(hierarchy.getPointer());

   void *p_pflotran_sim = NULL;

   f_create_simulation_(&p_pflotran_sim);

   f_set_hierarchy_ptr_(&p_pflotran_sim, &p_samr_hierarchy);

   f_initialize_simulation_(&p_pflotran_sim);

#if 0
   void *p_pflowhierarchy=NULL;
   f_create_hierarchy_data_(&p_pflowhierarchy);

#if 0
   void *p_integrator = NULL;
   f_create_integrator_(&p_integrator);
#endif

   void *p_samr_hierarchy=(void *)(hierarchy.getPointer());
   gridparameters *params = new gridparameters;

   hier::BoxArray<NDIM> physicalDomain = hierarchy->getGridGeometry()->getPhysicalDomain();

   // assume one box representing the domain for now
   hier::Box< NDIM > physicalBox = physicalDomain.getBox(0);
   params->p_grid = p_pflowhierarchy;
#if 0
   params->p_timestep = p_integrator;
#else
   params->p_timestep = NULL;
#endif

   params->igeom = 1;
   params->nx = physicalBox.numberCells(0);
   params->ny = physicalBox.numberCells(1);
   params->nz = physicalBox.numberCells(2);
   params->npx = 1;
   params->npy = 1;
   params->npz = 1;
   params->nphase = 1;
   params->nlevels = hierarchy->getNumberLevels();
   params->usesamrai = PETSC_TRUE;
   params->p_samr_hierarchy = p_samr_hierarchy;

#if 0
   // Create the application
   PIMSApplicationParameters* application_parameters = new PIMSApplicationParameters();
   application_parameters->d_hierarchy = hierarchy;
   application_parameters->d_db = input_db->getDatabase("PIMS");
   application_parameters->d_pims_parameters = params;

   PIMSApplication* application  = new PIMSApplication( application_parameters );
   
   // Initialize x0
   application->setInitialConditions(0.0);
   
   // Create time integrator
   algs::TimeIntegratorParameters* integration_parameters = new algs::TimeIntegratorParameters();
   tbox::Pointer< tbox::Database > integrator_db = input_db->getDatabase("PredictorCorrector");
   integration_parameters->d_db = integrator_db;
   integration_parameters->dt_method = 1;
   integration_parameters->d_ic_vector = application->get_x();
   integration_parameters->d_application_strategy = application;
   integration_parameters->d_integrator_name = "PredictorCorrector";
   algs::ForwardEulerTimeIntegrator* integrator = new algs::ForwardEulerTimeIntegrator(integration_parameters);

   // Create x(t)
   tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > x_t;
   x_t = integrator->getCurrentSolution();

   // Loop through time
   double dt = integrator->getCurrentDt();
   double time = integrator->getCurrentTime();
   tbox::pout << "dt = " << dt << "\n";
   int i = 1;
   int retval;
   double last_save = 0.0;
   char buffer [50];

   while ( time < integrator->getFinalTime() ) 
   {
      // Advance the solution
      if (i==1)
         integrator->advanceSolution(dt,true);
      else
         integrator->advanceSolution(dt,false);
      // Get the current solution and time
      x_t = integrator->getCurrentSolution();
      time = integrator->getCurrentTime();
      // Update the timestep
      dt = integrator->getNextDt(true,retval);
      // Write data
      if ( i%5==0 ) {
         sprintf(buffer,"t = %8.4f, dt = %8.6f\n",time,dt);
         tbox::pout << buffer;
      }

      i++;
   }

   // That's all, folks!
   delete application;
#endif
#endif

   /* 
    * That's all, folks!
    */
   PetscFinalize();
   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();

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

