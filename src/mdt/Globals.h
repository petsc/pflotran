#ifndef GLOBALS_H_
#define GLOBALS_H_

#define EAST 2 
#define WEST 4 
#define NORTH 3 
#define SOUTH 1 
#define TOP 6 
#define BOTTOM 5 

extern int myrank;
extern int commsize;

const double density0 = 998.32;
const double viscosity0 = 8.9e-4;
//const double gravity = 9.8068;
const double gravity = 9.81;
const double Ss = 1.e-6;  // specific storage (m^-1)
//const double gravity = 0.;
const double patm = 101325.;
const double p_threshhold = 0.03*gravity*density0;
const double permeability0 = 1.e-12;
const double porosity0 = 0.33;
const double betap = -4./(density0*gravity);

const double perm_to_hc = density0*gravity/viscosity0;
const double hc_to_perm = 1./perm_to_hc;
const double day_to_sec = 24.*3600.;

#define NULL_FLAG 0             ((unsigned long int) 0)        
#define WEST_DIR_WEST_FACE      ((unsigned long int) (1 << 0)) 
#define WEST_DIR_EAST_FACE      ((unsigned long int) (1 << 1)) 
#define WEST_DIR_SOUTH_FACE     ((unsigned long int) (1 << 2)) 
#define WEST_DIR_NORTH_FACE     ((unsigned long int) (1 << 3)) 
#define WEST_DIR_BOTTOM_FACE    ((unsigned long int) (1 << 4)) 
#define WEST_DIR_TOP_FACE       ((unsigned long int) (1 << 5)) 
#define EAST_DIR_WEST_FACE      ((unsigned long int) (1 << 6)) 
#define EAST_DIR_EAST_FACE      ((unsigned long int) (1 << 7)) 
#define EAST_DIR_SOUTH_FACE     ((unsigned long int) (1 << 8)) 
#define EAST_DIR_NORTH_FACE     ((unsigned long int) (1 << 9)) 
#define EAST_DIR_BOTTOM_FACE    ((unsigned long int) (1 << 10)) 
#define EAST_DIR_TOP_FACE       ((unsigned long int) (1 << 11)) 
#define SOUTH_DIR_WEST_FACE     ((unsigned long int) (1 << 12)) 
#define SOUTH_DIR_EAST_FACE     ((unsigned long int) (1 << 13)) 
#define SOUTH_DIR_SOUTH_FACE    ((unsigned long int) (1 << 14)) 
#define SOUTH_DIR_NORTH_FACE    ((unsigned long int) (1 << 15)) 
#define SOUTH_DIR_BOTTOM_FACE   ((unsigned long int) (1 << 16)) 
#define SOUTH_DIR_TOP_FACE      ((unsigned long int) (1 << 17)) 
#define NORTH_DIR_WEST_FACE     ((unsigned long int) (1 << 18)) 
#define NORTH_DIR_EAST_FACE     ((unsigned long int) (1 << 19)) 
#define NORTH_DIR_SOUTH_FACE    ((unsigned long int) (1 << 20)) 
#define NORTH_DIR_NORTH_FACE    ((unsigned long int) (1 << 21)) 
#define NORTH_DIR_BOTTOM_FACE   ((unsigned long int) (1 << 22)) 
#define NORTH_DIR_TOP_FACE      ((unsigned long int) (1 << 23)) 
#define BOTTOM_DIR_WEST_FACE    ((unsigned long int) (1 << 24)) 
#define BOTTOM_DIR_EAST_FACE    ((unsigned long int) (1 << 25)) 
#define BOTTOM_DIR_SOUTH_FACE   ((unsigned long int) (1 << 26)) 
#define BOTTOM_DIR_NORTH_FACE   ((unsigned long int) (1 << 27)) 
#define BOTTOM_DIR_BOTTOM_FACE  ((unsigned long int) (1 << 28)) 
#define BOTTOM_DIR_TOP_FACE     ((unsigned long int) (1 << 29)) 
#define TOP_DIR_TOP_FACE        ((unsigned long int) (1 << 30)) 
//#define TOP_DIR_WEST_FACE       ((unsigned long int) (1 << 31)) 
//#define TOP_DIR_EAST_FACE       ((unsigned long int) (1 << 32)) 
//#define TOP_DIR_SOUTH_FACE      ((unsigned long int) (1 << 33)) 
//#define TOP_DIR_NORTH_FACE      ((unsigned long int) (1 << 34)) 
//#define TOP_DIR_BOTTOM_FACE     ((unsigned long int) (1 << 35)) 

#endif /*GLOBALS_H_*/
