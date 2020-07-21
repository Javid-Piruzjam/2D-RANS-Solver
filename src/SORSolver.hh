#ifndef SOR_SOLVER_HH
#define SOR_SOLVER_HH

#include "StaggeredGrid.hh"


class SORSolver
{
public:
   SORSolver ( ) {}
   // Constructor to manually create SORSolver
   SORSolver ( int itermax , real eps , real omg );

   // Constructor to create a SORSolver from a parsed configuration file
   SORSolver ( const FileReader & configuration );

   real omg() {return omg_;}

   // solve the pressure equation on the staggered grid
   bool solve(StaggeredGrid & grid);


   void residual_calculation(StaggeredGrid & grid);

   real resL2Norm(StaggeredGrid & grid);


   void initGridSetup1(StaggeredGrid & grid);

   void initGridSetup2(StaggeredGrid & grid);

   void normalizer(StaggeredGrid & grid);

private:
   int itermax_;
   real eps_;
   real omg_;
   int checkfrequency_;
};




#endif //SOR_SOLVER_HH
