#ifndef TYPES_HH
#define TYPES_HH


// This typedef makes it possible to switch between float and double accuracy
// please do not use "float" or "double" directly in your code, but use real instead
typedef double real;



// Enumeration of boundary conditions
typedef enum { NOSLIP, SLIP, INFLOW, OUTFLOW, PERIODIC } BCTYPE;

//for directions
typedef enum {NORTH, EAST, SOUTH, WEST} Direction;


#endif //TYPES_HH
