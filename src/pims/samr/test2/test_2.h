//
// $Id: test_gradient.h 1832 2005-08-23 21:10:00Z bphilip $
// $Revision: 1832 $
// $Date: 2005-08-23 15:10:00 -0600 (Tue, 23 Aug 2005) $
//

/**
 * \file test_gradient.h
 * 
 * This simple example illustrates the computation of a gradient on an
 * AMR hierarchy.  
 */

#include <string>

#include "tbox/Database.h"
#include "PatchHierarchy.h"
#include "tbox/Pointer.h"
#include "StandardTagAndInitStrategy.h"

using namespace SAMRAI;

/**
 * Extract name of input and log files from command line.
 * \param argc        
 *        Number of command line arguments.
 * \param argv        
 *        Pointers to command line arguments.
 * \param input_file        
 *        Name of input file as it appears on command line.
 * \param log_file        
 *        Name of log file as it appears on command line.
 */
void processCommandLine(int argc,
                        char* argv[],
                        string& input_file,
                        string& log_file);

/**
 * Generate a patch hierarchy from information specified in input
 * database.
 * \param input_db        
 *        Input data base.
 * \param user_tagging_strategy        
 *        tbox::Pointer to user-provided cell tagging class.
 * \param hierarchy                    
 *        hier::Patch<NDIM> hierarchy initialized by this routine.
 */
void initializeAMRHierarchy(tbox::Pointer<tbox::Database> &input_db, 
			    mesh::StandardTagAndInitStrategy<NDIM>* user_tagging_strategy, 
			    tbox::Pointer<hier::PatchHierarchy<NDIM> > &hierarchy); 

/**
 * Allocate data for variables on a patch hierarchy.
 * \param hierarchy        
 *        hier::Patch<NDIM> hierarchy.
 * \param phi_id        
 *        Descriptor index of cell centered variable.
 * \param grad_phi_id        
 *        Descriptor index of face centered variable.
 */
void createVariableOnHierarchy(tbox::Pointer<hier::PatchHierarchy<NDIM> > &hierarchy,
                                tbox::Pointer< pdat::CellVariable<NDIM,double> > &phi,
				int &phi_id);

void
initializeVariableOnHierarchy(tbox::Pointer<hier::PatchHierarchy<NDIM> > &hierarchy,
                              const int phi_id);
