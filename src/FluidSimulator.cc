#include "FluidSimulator.hh"
#include "Debug.hh"
#include <math.h>
#include <algorithm>


FluidSimulator::FluidSimulator( const FileReader & conf ):grid_(StaggeredGrid(conf)), solver_(SORSolver(conf))
{
    gamma_ = conf.getRealParameter("gamma");
    Re_ = conf.getRealParameter("Re");
    gx_ = conf.getRealParameter("GX");
    gy_ = conf.getRealParameter("GY");
    dt_ = conf.getRealParameter("dt");
    U_INIT_ = conf.getRealParameter("U_INIT");
    V_INIT_ = conf.getRealParameter("V_INIT");
    P_INIT_ = conf.getRealParameter("P_INIT");
    K_INIT_ = conf.getRealParameter("K_INIT");
    E_INIT_ = conf.getRealParameter("E_INIT");

    c_e_ = conf.getRealParameter("C_e");
    c_mu_ = conf.getRealParameter("C_mu");
    c1_ = conf.getRealParameter("C1");
    c2_ = conf.getRealParameter("C2");

    safetyfactor_= conf.getRealParameter("safetyfactor");
    normalizationfrequency_ = conf.getIntParameter ("normalizationfrequency");
    outputinterval_= conf.getIntParameter ("outputinterval");
    
    real b = InletLength();
    nue_ = b * conf.getRealParameter("boundary_velocity_W")/Re_;   


    CHECK_MSG((gamma_>=0 && gamma_<=1), "Check the value of gamma!\n");

//defining boundary conditions and their value if provided!
    NB_=conf.getStringParameter("boundary_condition_N");
    NV_=conf.getRealParameter("boundary_velocity_N");
    SB_=conf.getStringParameter("boundary_condition_S");
    SV_=conf.getRealParameter("boundary_velocity_S");
    EB_=conf.getStringParameter("boundary_condition_E");
    EV_=conf.getRealParameter("boundary_velocity_E");
    WB_=conf.getStringParameter("boundary_condition_W");
    WV_=conf.getRealParameter("boundary_velocity_W");

}


real FluidSimulator::InletLength()
    {
      int n=0;
      for(int j=0 ; j<=grid_.jmax() ; j++)
	if(grid_.isFluid(1,j))
	  n++;
	
      return (n * grid_.dy());
    }

    

void FluidSimulator::simulate (real duration)
{

  unsigned int n=0;
  real t=0;
  VTKWriter vtkWriter (grid_, "lidDrivenCavity", true, true, false, false );
  
  //********************* Initialization **********************
  for(int j=0 ; j<=grid_.jmax() ; j++)
    if(grid_.isFluid(1,j))
      for (int i=0; i<=grid_.imax(); i++)
        grid_.u()(i,j)=1;
    else
      for (int i=0; i<=grid_.imax(); i++)
        grid_.u()(i,j)=0;

  grid_.valueInitialization(grid_.u(), U_INIT_);
  grid_.valueInitialization(grid_.v(), V_INIT_);
  grid_.valueInitialization(grid_.p(), P_INIT_);

  while (t<duration)
  {
    determineNextDT();
    refreshBoundaries();
    computeFG();
    composeRHS();
    solver_.solve(grid_);
    if (n%normalizationfrequency_==0)  solver_.normalizer(grid_);
    updateVelocities();
    if (n%outputinterval_==0)  vtkWriter.write();
    t += dt_;
    n += 1;
    if (n%10==0) std::cout<<"time step: "<<n<<"  time: "<<t<<'\n';
  }
}

void FluidSimulator::simulateTimeStepCount( unsigned int nrOfTimeSteps )
{

  unsigned int n=0;
  real t=0;
  VTKWriter vtkWriter (grid_, "lidDrivenCavity", true, true, false, false );

  //********************* Initialization **********************
  for(int j=0 ; j<=grid_.jmax() ; j++)
    if(grid_.isFluid(1,j))
      for (int i=0; i<=grid_.imax(); i++)
        grid_.u()(i,j)=1;
    else
      for (int i=0; i<=grid_.imax(); i++)
        grid_.u()(i,j)=0;

  grid_.valueInitialization(grid_.u(), U_INIT_);
  grid_.valueInitialization(grid_.v(), V_INIT_);
  grid_.valueInitialization(grid_.p(), P_INIT_);

  while (n<nrOfTimeSteps)
  {
    determineNextDT();
    refreshBoundaries(); 
    computeFG();  
    composeRHS();
    solver_.solve(grid_);
    if (n%normalizationfrequency_==0)  solver_.normalizer(grid_);
    updateVelocities();
    if (n%outputinterval_==0)  vtkWriter.write();
    t += dt_;
    n += 1;
    if (n%10==0) std::cout<<"time step: "<<n<<"  time: "<<t<<'\n';
  }
}


void FluidSimulator::simulateTurbulence ( unsigned int nrOfTimeSteps )
{

  std::cout<<"\n--------------------------- Turbulence Modeling ------------------------------- \n\n";
  unsigned int n=0;
  real t=0;
  VTKWriter vtkWriter (grid_, "TurbulenceBackstep", true, true, true, true );

  //********************* Initialization **********************
  for(int j=0 ; j<=grid_.jmax() ; j++)
    if(grid_.isFluid(1,j))
      for (int i=0; i<=grid_.imax(); i++)
        grid_.u()(i,j)=1;
    else
      for (int i=0; i<=grid_.imax(); i++)
        grid_.u()(i,j)=0;

  grid_.valueInitialization(grid_.u(), U_INIT_);
  grid_.valueInitialization(grid_.v(), V_INIT_);
  grid_.valueInitialization(grid_.p(), P_INIT_);
  grid_.valueInitialization(grid_.K(), K_INIT_);
  grid_.valueInitialization(grid_.E(), E_INIT_);
  
  
  computeDistanceFromWalls(); //Calculates each cell's distance from the nearest wall, including obstacle surfaces

  while (n<nrOfTimeSteps)
  {
    //determineNextDT();
    refreshBoundaries();
    computeRANS_Coefs();
    computeKE();  //for turbulence case
    computeTurbulentFG(); //for turbulence case
    composeRHS();
    solver_.solve(grid_);
    if (n%normalizationfrequency_==0)  solver_.normalizer(grid_);
    updateVelocities();
    if (n%outputinterval_==0)  vtkWriter.write();
    t += dt_;
    n += 1;
    if (n%10==0) std::cout<<"time step: "<<n<<"  time: "<<t<<'\n';
  }
}


void FluidSimulator::computeDistanceFromWalls()
{
   real r;
   grid_.d().fill(100);

   for(int j=1 ; j<=grid_.jmax() ;j++)
        for(int i=1; i<=grid_.imax() ; i++){
            r = (j-0-0.5)*grid_.dy();
            if (r<grid_.d()(i,j))  grid_.d()(i,j)=r;

            r = (grid_.jmax()-j+0.5)*grid_.dy();
            if (r<grid_.d()(i,j))  grid_.d()(i,j)=r;
        }

   for(int jo=1 ; jo<=grid_.jmax() ;jo++)
    for(int io=1; io<=grid_.imax() ; io++)
      if(!grid_.isFluid(io,jo)){
	grid_.d()(io,jo)=0;
          for(int j=1 ; j<=grid_.jmax() ;j++){
            for(int i=1; i<=grid_.imax() ; i++){
		r = sqrt( ((i-io)*grid_.dx())*((i-io)*grid_.dx()) + ((j-jo)*grid_.dy())*((j-jo)*grid_.dy()) );
                if (r<grid_.d()(i,j))  grid_.d()(i,j)=r;
	    }
	  }
       }
}

void FluidSimulator::computeFG()
{
real Du2Dx;
real DuvDy;
real D2uDx2;
real D2uDy2;
//real DpDx;

real DuvDx;
real Dv2Dy;
real D2vDx2;
real D2vDy2;
//real DpDy;

for (int i=1 ; i<=grid_.imax() ; i++ )
    for (int j=1 ; j<=grid_.jmax() ; j++ )

    if(grid_.isFluid(i,j))
    {
      if (i != grid_.imax() ) {

        Du2Dx = (0.25/grid_.dx()) * ( (grid_.u()(i,j) + grid_.u(i+1,j,EAST))*(grid_.u()(i,j) + grid_.u(i+1,j,EAST))
                                    - (grid_.u(i-1,j,WEST) + grid_.u()(i,j))*(grid_.u(i-1,j,WEST) + grid_.u()(i,j)))
            + (0.25*gamma_/grid_.dx()) * ( fabs(grid_.u()(i,j) + grid_.u(i+1,j,EAST))*(grid_.u()(i,j) - grid_.u(i+1,j,EAST))
                                    - fabs(grid_.u(i-1,j,WEST) + grid_.u()(i,j))*(grid_.u(i-1,j,WEST) - grid_.u()(i,j)));


        DuvDy = (0.25/grid_.dy()) * ( (grid_.v()(i,j) + grid_.v(i+1,j,EAST))*(grid_.u()(i,j) + grid_.u(i,j+1,NORTH))
                                    - (grid_.v(i,j-1,SOUTH) + grid_.v(i+1,j-1,SOUTH))*(grid_.u(i,j-1,SOUTH) + grid_.u()(i,j)))
            + (0.25*gamma_/grid_.dy()) * ( fabs(grid_.v()(i,j) + grid_.v(i+1,j,EAST))*(grid_.u()(i,j) - grid_.u(i,j+1,NORTH))
                                    - fabs(grid_.v(i,j-1,SOUTH) + grid_.v(i+1,j-1,SOUTH))*(grid_.u(i,j-1,SOUTH) - grid_.u()(i,j)));

        D2uDx2 = (1./((grid_.dx())*(grid_.dx()))) * ( grid_.u(i+1,j,EAST) - 2*grid_.u()(i,j) + grid_.u(i-1,j,WEST) );
        D2uDy2 = (1./((grid_.dy())*(grid_.dy()))) * ( grid_.u(i,j+1,NORTH) - 2*grid_.u()(i,j) + grid_.u(i,j-1,SOUTH) );
        //DpDx = (1/(grid_.dx())) * ( grid_.p()(i+1,j) - grid_.p()(i,j) );

        grid_.f()(i,j) = grid_.u()(i,j) + dt_ * ( (1./Re_)*(D2uDx2 + D2uDy2) - Du2Dx - DuvDy + gx_  );

        }


      if(j != grid_.jmax() ) {

        DuvDx = (0.25/grid_.dx()) * ( (grid_.u()(i,j) + grid_.u(i,j+1,NORTH))*(grid_.v()(i,j) + grid_.v(i+1,j,EAST))
                                    - (grid_.u(i-1,j,WEST) + grid_.u(i-1,j+1,NORTH))*(grid_.v(i-1,j,WEST) + grid_.v()(i,j)))
            + (0.25*gamma_/grid_.dx()) * (fabs(grid_.u()(i,j) + grid_.u(i,j+1,NORTH))*(grid_.v()(i,j) - grid_.v(i+1,j,EAST))
                                    - fabs(grid_.u(i-1,j,WEST) + grid_.u(i-1,j+1,NORTH))*(grid_.v(i-1,j,WEST) - grid_.v()(i,j)));

        Dv2Dy = (0.25/grid_.dy()) * ( (grid_.v()(i,j) + grid_.v(i,j+1,NORTH))*(grid_.v()(i,j) + grid_.v(i,j+1,NORTH))
                                    - (grid_.v(i,j-1,SOUTH) + grid_.v()(i,j))*(grid_.v(i,j-1,SOUTH) + grid_.v()(i,j)))
            + (0.25*gamma_/grid_.dy()) * ( (fabs(grid_.v()(i,j) + grid_.v(i,j+1,NORTH)))*(grid_.v()(i,j) - grid_.v(i,j+1,NORTH))
                                    - fabs(grid_.v(i,j-1,SOUTH) + grid_.v()(i,j))*(grid_.v(i,j-1,SOUTH) - grid_.v()(i,j)));

        D2vDx2 = (1/((grid_.dx())*(grid_.dx()))) * ( grid_.v(i+1,j,EAST) - 2*grid_.v()(i,j) + grid_.v(i-1,j,WEST) );
        D2vDy2 = (1/((grid_.dy())*(grid_.dy()))) * ( grid_.v(i,j+1,NORTH) - 2*grid_.v()(i,j) + grid_.v(i,j-1,SOUTH) );
        //DpDy = (1/(grid_.dy())) * ( grid_.p()(i,j+1) - grid_.p()(i,j) );

        grid_.g()(i,j) = grid_.v()(i,j) + dt_ * (  (1/Re_)*(D2vDx2 + D2vDy2) - DuvDx - Dv2Dy + gy_   );
      }
    }


    for (int i=1 ; i<=grid_.imax() ; i++) {

        grid_.g()(i,0) = grid_.v()(i,0);
        grid_.g()(i,grid_.jmax()) = grid_.v()(i,grid_.jmax());

    }
    for (int j=1 ; j<=grid_.jmax() ; j++){
        grid_.f()(0,j) = grid_.u()(0,j);
        grid_.f()(grid_.imax(),j) = grid_.u()(grid_.imax(),j);
    }
}

void FluidSimulator::computeTurbulentFG()
{
real termF1, termF2, DkDx, Du2Dx, DuvDy;
real termG1, termG2, DkDy, DuvDx, Dv2Dy;

for (int i=1 ; i<=grid_.imax() ; i++ )
    for (int j=1 ; j<=grid_.jmax() ; j++ )

    if(grid_.isFluid(i,j))
    {
      if (i != grid_.imax() ) {

        termF1 = (1./(grid_.dx()*grid_.dx())) * ( (nue_+grid_.nueT()(i+1,j))*(grid_.u(i+1,j,EAST)-grid_.u()(i,j))
                                                - (nue_+grid_.nueT()(i,j))*(grid_.u()(i,j)-grid_.u(i-1,j,WEST)) );

        termF2 =(1./(grid_.dy())) * ((nue_ + 0.25*(grid_.nueT()(i+1,j+1)+grid_.nueT()(i,j+1)+grid_.nueT()(i+1,j)+grid_.nueT()(i,j)) )
                                   *( ( grid_.u(i,j+1,NORTH)-grid_.u()(i,j) )/grid_.dy() + ( grid_.v(i+1,j,EAST)-grid_.v()(i,j) )/grid_.dx() )
                                   - (nue_ + 0.25*(grid_.nueT()(i+1,j)+grid_.nueT()(i,j)+grid_.nueT()(i+1,j-1)+grid_.nueT()(i,j-1)) )
                                   *( ( grid_.u()(i,j)-grid_.u(i,j-1,SOUTH) )/grid_.dy() + ( grid_.v(i+1,j-1,EAST)-grid_.v(i,j-1,SOUTH) )/grid_.dx() ));

        DkDx = (1./(grid_.dx())) * ( grid_.K(i+1,j,EAST) - grid_.K()(i,j) );

        Du2Dx = (0.25/grid_.dx()) * ( (grid_.u()(i,j) + grid_.u(i+1,j,EAST))*(grid_.u()(i,j) + grid_.u(i+1,j,EAST))
                                    - (grid_.u(i-1,j,WEST) + grid_.u()(i,j))*(grid_.u(i-1,j,WEST) + grid_.u()(i,j)))
                + (0.25*gamma_/grid_.dx()) * ( fabs(grid_.u()(i,j) + grid_.u(i+1,j,EAST))*(grid_.u()(i,j) - grid_.u(i+1,j,EAST))
                                    - fabs(grid_.u(i-1,j,WEST) + grid_.u()(i,j))*(grid_.u(i-1,j,WEST) - grid_.u()(i,j)));

        DuvDy = (0.25/grid_.dy()) * ( (grid_.v()(i,j) + grid_.v(i+1,j,EAST))*(grid_.u()(i,j) + grid_.u(i,j+1,NORTH))
                                    - (grid_.v(i,j-1,SOUTH) + grid_.v(i+1,j-1,SOUTH))*(grid_.u(i,j-1,SOUTH) + grid_.u()(i,j)))
            + (0.25*gamma_/grid_.dy()) * ( fabs(grid_.v()(i,j) + grid_.v(i+1,j,EAST))*(grid_.u()(i,j) - grid_.u(i,j+1,NORTH))
                                    - fabs(grid_.v(i,j-1,SOUTH) + grid_.v(i+1,j-1,SOUTH))*(grid_.u(i,j-1,SOUTH) - grid_.u()(i,j)));

        grid_.f()(i,j) = grid_.u()(i,j) + dt_ * (2*termF1 + termF2 - (2.0/3.0)*DkDx - Du2Dx - DuvDy  + gx_  );
        }


      if(j != grid_.jmax() ) {

        termG1 =(1./(grid_.dx())) * ((nue_ + 0.25*(grid_.nueT()(i+1,j+1)+grid_.nueT()(i,j+1)+grid_.nueT()(i+1,j)+grid_.nueT()(i,j)) )
                                   *( ( grid_.v(i+1,j,EAST)-grid_.v()(i,j) )/grid_.dx() + ( grid_.u(i,j+1,NORTH)-grid_.u()(i,j) )/grid_.dy() )
                                   - (nue_ + 0.25*(grid_.nueT()(i,j+1)+grid_.nueT()(i-1,j+1)+grid_.nueT()(i,j)+grid_.nueT()(i-1,j)) )
                                   *( ( grid_.v()(i,j)-grid_.v(i-1,j,WEST) )/grid_.dx() + ( grid_.u(i-1,j+1,NORTH)-grid_.u(i-1,j,WEST) )/grid_.dy() ) );

        termG2 =  (1./(grid_.dy()*grid_.dy())) * ( (nue_+grid_.nueT()(i,j+1))*(grid_.v(i,j+1,NORTH)-grid_.v()(i,j))
                                                - (nue_+grid_.nueT()(i,j))*(grid_.v()(i,j)-grid_.v(i,j-1,SOUTH)) );

        DkDy = (1./(grid_.dy())) * ( grid_.K(i,j+1,NORTH) - grid_.K()(i,j) );

        DuvDx = (0.25/grid_.dx()) * ( (grid_.u()(i,j) + grid_.u(i,j+1,NORTH))*(grid_.v()(i,j) + grid_.v(i+1,j,EAST))
                                    - (grid_.u(i-1,j,WEST) + grid_.u(i-1,j+1,NORTH))*(grid_.v(i-1,j,WEST) + grid_.v()(i,j)))
            + (0.25*gamma_/grid_.dx()) * (fabs(grid_.u()(i,j) + grid_.u(i,j+1,NORTH))*(grid_.v()(i,j) - grid_.v(i+1,j,EAST))
                                    - fabs(grid_.u(i-1,j,WEST) + grid_.u(i-1,j+1,NORTH))*(grid_.v(i-1,j,WEST) - grid_.v()(i,j)));

        Dv2Dy = (0.25/grid_.dy()) * ( (grid_.v()(i,j) + grid_.v(i,j+1,NORTH))*(grid_.v()(i,j) + grid_.v(i,j+1,NORTH))
                                    - (grid_.v(i,j-1,SOUTH) + grid_.v()(i,j))*(grid_.v(i,j-1,SOUTH) + grid_.v()(i,j)))
            + (0.25*gamma_/grid_.dy()) * ( (fabs(grid_.v()(i,j) + grid_.v(i,j+1,NORTH)))*(grid_.v()(i,j) - grid_.v(i,j+1,NORTH))
                                    - fabs(grid_.v(i,j-1,SOUTH) + grid_.v()(i,j))*(grid_.v(i,j-1,SOUTH) - grid_.v()(i,j)));

        grid_.g()(i,j) = grid_.v()(i,j) + dt_ * ( termG1 + 2*termG2 - (2.0/3.0)*DkDy - DuvDx - Dv2Dy + gy_ );
      }
    }


    for (int i=1 ; i<=grid_.imax() ; i++) {

        grid_.g()(i,0) = grid_.v()(i,0);
        grid_.g()(i,grid_.jmax()) = grid_.v()(i,grid_.jmax());

    }
    for (int j=1 ; j<=grid_.jmax() ; j++){
        grid_.f()(0,j) = grid_.u()(0,j);
        grid_.f()(grid_.imax(),j) = grid_.u()(grid_.imax(),j);
    }
}


void FluidSimulator::computeRANS_Coefs()
{
 real Re_t, Re_d;
 
 for (int i=1 ; i<=grid_.imax() ; i++ )
    for (int j=1 ; j<=grid_.jmax() ; j++ )
       if(grid_.isFluid(i,j))
       {
        Re_t =  grid_.K()(i,j)*grid_.K()(i,j) / (nue_ * grid_.E()(i,j));
        Re_d = sqrt(grid_.K()(i,j)) * grid_.d()(i,j) / nue_;
        grid_.f_mu()(i,j) = (1-exp(-0.0165*Re_d))*(1-exp(-0.0165*Re_d)) * (1 + 20.5/Re_t);

        grid_.f1()(i,j) = 1 + (0.05/grid_.f_mu()(i,j))*(0.05/grid_.f_mu()(i,j))*(0.05/grid_.f_mu()(i,j));
	
        grid_.f2()(i,j) = 1 - exp(-Re_t*Re_t);

        grid_.nueT()(i,j) = c_mu_ * grid_.f_mu()(i,j) * grid_.K()(i,j)*grid_.K()(i,j) / grid_.E()(i,j);
       }

}

void FluidSimulator::computeKE()
{

for (int i=1 ; i<=grid_.imax() ; i++ )
    for (int j=1 ; j<=grid_.jmax() ; j++ )

    if(grid_.isFluid(i,j))
    {
    // Calculations for K:
    real termK1 = (1./(grid_.dx()*grid_.dx())) * (0.5*(grid_.nueT()(i,j)+grid_.nueT()(i+1,j))*(grid_.K(i+1,j,EAST)-grid_.K()(i,j))
                                                - 0.5*(grid_.nueT()(i-1,j)+grid_.nueT()(i,j))*(grid_.K()(i,j)-grid_.K(i-1,j,WEST)));
    real termK2 = (1./(grid_.dy()*grid_.dy())) * (0.5*(grid_.nueT()(i,j)+grid_.nueT()(i,j+1))*(grid_.K(i,j+1,NORTH)-grid_.K()(i,j))
                                                - 0.5*(grid_.nueT()(i,j-1)+grid_.nueT()(i,j))*(grid_.K()(i,j)-grid_.K(i,j-1,SOUTH)));
    real DukDx = (0.5/grid_.dx())*(grid_.u()(i,j)*(grid_.K()(i,j)+grid_.K(i+1,j,EAST)) - (grid_.u(i-1,j,WEST)*(grid_.K(i-1,j,WEST)+grid_.K()(i,j))))
                  +((0.5*gamma_)/grid_.dx())*(fabs(grid_.u()(i,j))*(grid_.K()(i,j)-grid_.K(i+1,j,EAST)) - (fabs(grid_.u(i-1,j,WEST))*(grid_.K(i-1,j,WEST)-grid_.K()(i,j))));
    real DvkDy =(0.5/grid_.dy())*(grid_.v()(i,j)*(grid_.K()(i,j)+grid_.K(i,j+1,NORTH)) - (grid_.v(i,j-1,SOUTH)*(grid_.K(i,j-1,SOUTH)+grid_.K()(i,j))))
                  +((0.5*gamma_)/grid_.dy())*(fabs(grid_.v()(i,j))*(grid_.K()(i,j)-grid_.K(i,j+1,NORTH)) - (fabs(grid_.v(i,j-1,SOUTH))*(grid_.K(i,j-1,SOUTH)-grid_.K()(i,j))));
    real DuDx = (1./grid_.dx())*(grid_.u()(i,j)-grid_.u(i-1,j,WEST));
    real DuDy = (0.25/grid_.dy())*(grid_.u(i,j+1,NORTH)+grid_.u(i-1,j+1,NORTH)-grid_.u(i,j-1,SOUTH)-grid_.u(i-1,j-1,SOUTH));
    real DvDx = (0.25/grid_.dx())*(grid_.v(i+1,j,EAST)+grid_.v(i+1,j-1,SOUTH)-grid_.v(i-1,j,WEST)-grid_.v(i-1,j-1,SOUTH));
    real DvDy  = (1./grid_.dy())*(grid_.v()(i,j)-grid_.v(i,j-1,SOUTH));
    real grads =4*(DuDx*DuDx) + 2*(DuDy+DvDx)*(DuDy+DvDx)+ 4*(DvDy*DvDy);
    real termK3 = 0.5 * grid_.nueT()(i,j) * grads;
    real K_Old= grid_.K()(i,j);
    grid_.K()(i,j) = grid_.K()(i,j) + dt_ * (termK1 + termK2 - DukDx - DvkDy + termK3 - grid_.E()(i,j));


    // Calculations for Epsilon:
    real termE1 = (1./(grid_.dx()*grid_.dx())) * (0.25*(grid_.nueT()(i,j)+grid_.nueT()(i+1,j))*(grid_.f_mu()(i,j)+grid_.f_mu()(i+1,j))*(grid_.E(i+1,j,EAST)-grid_.E()(i,j))
                                               - 0.25*(grid_.nueT()(i-1,j)+grid_.nueT()(i,j))*(grid_.f_mu()(i-1,j)+grid_.f_mu()(i,j))*(grid_.E()(i,j)-grid_.E(i-1,j,WEST)));
    real termE2 = (1./(grid_.dy()*grid_.dy())) * (0.25*(grid_.nueT()(i,j)+grid_.nueT()(i,j+1))*(grid_.f_mu()(i,j)+grid_.f_mu()(i,j+1))*(grid_.E(i,j+1,NORTH)-grid_.E()(i,j))
                                             - 0.25*(grid_.nueT()(i,j-1)+grid_.nueT()(i,j))*(grid_.f_mu()(i,j-1)+grid_.f_mu()(i,j))*(grid_.E()(i,j)-grid_.E(i,j-1,SOUTH)));

    real DueDx = (0.5/grid_.dx())*(grid_.u()(i,j)*(grid_.E()(i,j)+grid_.E(i+1,j,EAST)) - (grid_.u(i-1,j,WEST)*(grid_.E(i-1,j,WEST)+grid_.E()(i,j))))
                  +((0.5*gamma_)/grid_.dx())*(fabs(grid_.u()(i,j))*(grid_.E()(i,j)-grid_.E(i+1,j,EAST)) - (fabs(grid_.u(i-1,j,WEST))*(grid_.E(i-1,j,WEST)-grid_.E()(i,j))));
    real DveDy =(0.5/grid_.dy())*(grid_.v()(i,j)*(grid_.E()(i,j)+grid_.E(i,j+1,NORTH)) - (grid_.v(i,j-1,SOUTH)*(grid_.E(i,j-1,SOUTH)+grid_.E()(i,j))))
                  +((0.5*gamma_)/grid_.dy())*(fabs(grid_.v()(i,j))*(grid_.E()(i,j)-grid_.E(i,j+1,NORTH)) - (fabs(grid_.v(i,j-1,SOUTH))*(grid_.E(i,j-1,SOUTH)-grid_.E()(i,j))));
		  
    real termE3 = 0.5 * c1_ * grid_.f1()(i,j)* K_Old* grads;
    real termE4 = c2_ * grid_.f2()(i,j) * grid_.E()(i,j)*grid_.E()(i,j) / K_Old;
    grid_.E()(i,j) = grid_.E()(i,j) + dt_ * ((c_e_/c_mu_)*(termE1 + termE2) - DueDx - DveDy + termE3 - termE4);
    }
}


void FluidSimulator::composeRHS()
{
    for (int i=1 ; i<=grid_.imax() ; i++ )
        for (int j=1 ; j<=grid_.jmax() ; j++ )
          if(grid_.isFluid(i,j))
            grid_.rhs()(i,j)= (1/dt_) * ((grid_.f()(i,j)-grid_.f()(i-1,j))/grid_.dx() + (grid_.g()(i,j)-grid_.g()(i,j-1))/grid_.dy());
}

void FluidSimulator::updateVelocities()
{
    for (int i=1 ; i<grid_.imax() ; i++ )
        for (int j=1 ; j<=grid_.jmax() ; j++ )
            if(grid_.isFluid(i,j)){
               if(grid_.isFluid(i+1,j))
                 grid_.u()(i,j) = grid_.f()(i,j) - (dt_/grid_.dx())*(grid_.p(i+1,j,EAST)-grid_.p()(i,j));
               else
                grid_.u()(i,j) = 0;
            }

    for (int i=1 ; i<=grid_.imax() ; i++ )
        for (int j=1 ; j<grid_.jmax() ; j++ )
            if(grid_.isFluid(i,j)){
               if(grid_.isFluid(i,j+1))
                 grid_.v()(i,j) = grid_.g()(i,j) - (dt_/grid_.dy())*(grid_.p(i,j+1,NORTH)-grid_.p()(i,j));
               else
                grid_.v()(i,j) = 0;
            }
}

void FluidSimulator::determineNextDT()
{
    if (safetyfactor_>0){
        real umax=0, vmax=0;
        real dx2 = grid_.dx() * grid_.dx();
        real dy2 = grid_.dy() * grid_.dy();
        real dt=1;

        for (int j=1 ; j<= grid_.jmax() ; j++)
            for (int i=1 ; i<= grid_.imax() ; i++)
                if(grid_.isFluid(i,j))
                umax = std::max (umax , grid_.u()(i,j));

        for (int j=1 ; j<= grid_.jmax() ; j++)
            for (int i=1 ; i<= grid_.imax() ; i++)
                if(grid_.isFluid(i,j))
                vmax = std::max (vmax , grid_.v()(i,j));

        if (umax!=0  &&  vmax!= 0 )   dt = std::min( (grid_.dx()/fabs(umax)),  (grid_.dy()/fabs(vmax)) );
        dt_ = safetyfactor_ * std::min( (0.5*Re_*((dx2*dy2)/(dx2+dy2))) , dt );
    }

}


BCTYPE FluidSimulator::BCtoType(const std::string BC)
{
  if (BC == "noslip") return NOSLIP;
  else if (BC == "NOSLIP") return NOSLIP;
  else if (BC == "outflow" ) return OUTFLOW;
  else if (BC == "inflow" ) return INFLOW;
  return NOSLIP;
}


void FluidSimulator::refreshBoundaries()
{
int imax = grid_.imax();
int jmax = grid_.jmax();

BCTYPE NorthBoundary = BCtoType(NB_);
BCTYPE SouthBoundary = BCtoType(SB_);
BCTYPE EastBoundary = BCtoType(EB_);
BCTYPE WestBoundary = BCtoType(WB_);

//NOSLIP=0, SLIP=1, INFLOW=2, OUTFLOW=3, PERIODIC=4

switch(NorthBoundary)
{
    case NOSLIP:
    {
        for(int i=0 ; i<=imax ; i++)
            grid_.u()(i,jmax+1)= 2*NV_ - grid_.u()(i,jmax);
        for(int i=1 ; i<=imax ; i++)
            grid_.v()(i,jmax)= 0;

        for(int i=0 ; i<=imax ; i++){
            grid_.K()(i,jmax+1) = NV_*NV_ -grid_.K()(i,jmax);
	    //grid_.K()(i,jmax+1) = grid_.K()(i,jmax);
            grid_.E()(i,jmax+1) = grid_.E()(i,jmax);
	    //grid_.E()(i,jmax+1) = -grid_.E()(i,jmax);
        }
    }
    break;

    case OUTFLOW:
    {
        for(int i=0 ; i<=imax ; i++)
            grid_.u()(i,jmax+1)= grid_.u()(i,jmax);
        for(int i=1 ; i<=imax ; i++)
            grid_.v()(i,jmax)= 0;

        for(int i=0 ; i<=imax ; i++){
            grid_.K()(i,jmax+1) = grid_.K()(i,jmax);
            grid_.E()(i,jmax+1) = grid_.E()(i,jmax);
        }
    }
    break;

    case INFLOW:
    {
        for(int i=0 ; i<=imax ; i++)
            grid_.u()(i,jmax+1)= - grid_.u()(i,jmax);
        for(int i=1 ; i<=imax ; i++)
            grid_.v()(i,jmax)= NV_;

        real K_in = 0.003 * NV_*NV_;
        real E_in = (c_mu_ * sqrt(K_in*K_in*K_in) ) / (0.03 * grid_.dx()*grid_.imax() ) ;
        for(int i=0 ; i<=imax ; i++){
            grid_.K()(i,jmax+1) = 2*K_in - grid_.K()(i,jmax);
            grid_.E()(i,jmax+1) = 2*E_in - grid_.E()(i,jmax);
        }
    }
    break;

    default:
    {
    std::cout<<"check the boundary condition type!\n";
    }
    break;
}


switch(SouthBoundary)
{
    case NOSLIP:
    {
        for(int i=0 ; i<=imax ; i++)
            grid_.u()(i,0)= 2*SV_ - grid_.u()(i,1);
        for(int i=1 ; i<=imax ; i++)
            grid_.v()(i,0)= 0;

        for(int i=0 ; i<=imax ; i++){
            grid_.K()(i,0) = SV_*SV_ - grid_.K()(i,1);
            grid_.E()(i,0) = grid_.E()(i,1);
        }
    }
    break;

    case OUTFLOW:
    {
       for(int i=0 ; i<=imax ; i++)
            grid_.u()(i,0)= grid_.u()(i,1);
       for(int i=1 ; i<=imax ; i++)
            grid_.v()(i,0)=grid_.v()(i,1);

       for(int i=0 ; i<=imax ; i++){
            grid_.K()(i,0) = grid_.K()(i,1);
            grid_.E()(i,0) = grid_.E()(i,1);
        }
    }
    break;

    case INFLOW:
    {
        for(int i=0 ; i<=imax ; i++)
            grid_.u()(i,0)= - grid_.u()(i,1);
        for(int i=1 ; i<=imax ; i++)
            grid_.v()(i,0)= SV_;

        real K_in = 0.003 * SV_*SV_;
        real E_in = (c_mu_ * sqrt(K_in*K_in*K_in) ) / (0.03 * grid_.dx()*grid_.imax() ) ;
        for(int i=0 ; i<=imax ; i++){
            grid_.K()(i,0) = 2*K_in - grid_.K()(i,1);
            grid_.E()(i,0) = 2*E_in - grid_.E()(i,1);
        }
    }
    break;

    default:
    {
    std::cout<<"check the boundary condition type!\n";
    }
    break;
}

switch(EastBoundary)
{
    case NOSLIP:
    {
        for(int j=0 ; j<=jmax ; j++)
          if (grid_.isFluid(imax,j))
            grid_.v()(imax+1,j)= 2*EV_ - grid_.v()(imax,j);
        for(int j=1 ; j<=jmax ; j++)
          if (grid_.isFluid(imax,j))
            grid_.u()(imax,j)= 0;

        for(int j=0 ; j<=jmax ; j++)
            if (grid_.isFluid(imax,j)){
            grid_.K()(imax+1,j) = EV_*EV_ - grid_.K()(imax,j);
            grid_.E()(imax+1,j) = grid_.E()(imax,j);
        }
    }
    break;

    case OUTFLOW:
    {
       for(int j=0 ; j<=jmax ; j++)
        if (grid_.isFluid(imax,j))
            grid_.v()(imax+1,j)= grid_.v()(imax,j);
       for(int j=1 ; j<=jmax ; j++)
        if (grid_.isFluid(imax,j))
            grid_.u()(imax,j)= grid_.u()(imax-1,j);

       for(int j=1 ; j<=jmax ; j++)
        if (grid_.isFluid(imax,j)){
            grid_.K()(imax+1,j) = grid_.K()(imax,j);
            grid_.E()(imax+1,j) = grid_.E()(imax,j);
        }
    }
    break;

    case INFLOW:
    {
       for(int j=0 ; j<=jmax ; j++)
        if (grid_.isFluid(imax,j))
            grid_.v()(imax+1,j)= - grid_.v()(imax,j);
       for(int j=1 ; j<=jmax ; j++)
        if (grid_.isFluid(imax,j))
            grid_.u()(imax,j)= EV_;


       real K_in = 0.003 * EV_*EV_;
       real E_in = (c_mu_ * sqrt(K_in*K_in*K_in) ) / (0.03 * grid_.dy()*grid_.jmax() ) ;
       for(int j=0 ; j<=jmax ; j++)
        if (grid_.isFluid(imax,j)){
            grid_.K()(imax+1,j) = 2*K_in - grid_.K()(imax,j);
            grid_.E()(imax+1,j) = 2*E_in - grid_.E()(imax,j);
        }

    }
    break;

    default:
    {
    std::cout<<"check the boundary condition type!\n";
    }
    break;
}


switch(WestBoundary)
{
    case NOSLIP:
    {
       for(int j=0 ; j<=jmax ; j++)
        if (grid_.isFluid(1,j))
            grid_.v()(0,j)= 2*WV_ - grid_.v()(1,j);
       for(int j=1 ; j<=jmax ; j++)
        if (grid_.isFluid(1,j))
            grid_.u()(0,j)= 0;

       for(int j=1 ; j<=jmax ; j++)
            if (grid_.isFluid(1,j)){
            grid_.K()(0,j) = WV_*WV_ - grid_.K()(1,j);
            grid_.E()(0,j) = grid_.E()(1,j);
        }
    }
    break;

    case OUTFLOW:
    {
       for(int j=0 ; j<=jmax ; j++)
        if (grid_.isFluid(1,j))
            grid_.v()(0,j)= grid_.v()(1,j);
       for(int j=1 ; j<=jmax ; j++)
        if (grid_.isFluid(1,j))
            grid_.u()(0,j)= grid_.u()(1,j);

       for(int j=0 ; j<=jmax ; j++)
        if (grid_.isFluid(1,j)){
            grid_.K()(0,j) = grid_.K()(1,j);
            grid_.E()(0,j) = grid_.E()(1,j);
        }
    }
    break;

    case INFLOW:
    {
      
       real b = InletLength();
       
      
       for(int j=0 ; j<=jmax ; j++)
         if (grid_.isFluid(1,j))
            grid_.v()(0,j)= - grid_.v()(1,j);
       for(int j=1 ; j<=jmax ; j++)
         if (grid_.isFluid(1,j))
            grid_.u()(0,j)= WV_;

       real K_in = 0.003 * WV_*WV_;
       real E_in = (c_mu_ * sqrt(K_in*K_in*K_in) ) / (0.03 * b ) ;
       for(int j=1 ; j<=jmax ; j++)
        if (grid_.isFluid(1,j)){
            grid_.K()(0,j) = 2*K_in - grid_.K()(1,j);
            grid_.E()(0,j) = 2*E_in - grid_.E()(1,j);
        }
    }
    break;

    default:
    {
    std::cout<<"check the boundary condition type!\n";
    }
    break;
  }
}
