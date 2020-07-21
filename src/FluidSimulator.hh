#ifndef __FLUID_SIMULATOR_H__
#define __FLUID_SIMULATOR_H__

#include "StaggeredGrid.hh"
#include "SORSolver.hh"
#include "VTKWriter.hh"

class FluidSimulator
{
  public:
      FluidSimulator( const FileReader & conf );

      /// Simulates a given time-length
      void simulate             ( real duration              );
      void simulateTimeStepCount( unsigned int nrOfTimeSteps );
      void simulateTurbulence( unsigned int nrOfTimeSteps );

      // Getter functions for the internally stored StaggeredGrid
            StaggeredGrid & grid() {return grid_;}
      const StaggeredGrid & grid() const {return grid_;}
            SORSolver & solver() {return solver_;}
      const SORSolver & solver() const {return solver_;}

      //void computeKE(U,V,KA,EP,FLAG,imax,jmax,delt,delx,dely,GX,GY,Re,gamma);
      real InletLength();
      void computeDistanceFromWalls();
      void computeRANS_Coefs();
      void computeKE();
      void composeRHS();
      void updateVelocities();
      void determineNextDT();
      BCTYPE BCtoType(const std::string);
      void refreshBoundaries();

  private:
      void computeFG();
      void computeTurbulentFG();

      StaggeredGrid grid_;
      SORSolver solver_;
      real gamma_;
      real Re_;
      real c_e_;
      real c_mu_;
      real c1_;
      real c2_;
      real nue_;
      real gx_;
      real gy_;
      real dt_;
      real Dt_;
      real U_INIT_;
      real V_INIT_;
      real P_INIT_;
      real K_INIT_;  //initial value of turbulent kinetic energy
      real E_INIT_;  //initial value of dissipation rate
      real safetyfactor_;
      int outputinterval_;
      int normalizationfrequency_;


      std::string NB_ ;
      std::string SB_ ;
      std::string EB_ ;
      std::string WB_;
      real NV_ = 0;
      real SV_ = 0;
      real EV_ = 0;
      real WV_ = 0;

};

#endif
