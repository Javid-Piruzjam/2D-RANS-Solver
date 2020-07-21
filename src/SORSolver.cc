#include "SORSolver.hh"
#include <iostream>
#include <math.h>
#include "Debug.hh"

#define PI 3.14159265

//class constructor using values
SORSolver::SORSolver ( int itermax , real eps , real omg )
{
    itermax_ = itermax;
    eps_ = eps;
    omg_ = omg;

    CHECK_MSG((omg_>=0 && omg_<=2), "Check the value of omaga!\n");
}
//class constructor using fileReader
SORSolver::SORSolver ( const FileReader & configuration )
{
    itermax_ = configuration.getIntParameter("itermax");
    eps_ = configuration.getRealParameter("eps");
    omg_ = configuration.getRealParameter("omg");
    checkfrequency_= configuration.getIntParameter("checkfrequency");

    CHECK_MSG((omg_>=0 && omg_<=2), "Check the value of omaga!\n");
}

bool SORSolver::solve( StaggeredGrid & grid )
{
  int iter = 0;
  real residual = eps_ + 1;
  //int imax = grid.imax();
  //int jmax = grid.imax();
  int imax = grid.imax();
  int jmax = grid.jmax();
  real dx2 = grid.dx()*grid.dx();
  real dy2 = grid.dy()*grid.dy();
  real coef = omg_* 0.5*((dx2*dy2)/(dx2+dy2));

  //******************* initialization method ************************
    //grid.randomInitialization();
  //  grid.randomInitialization2(grid.p());  // initilize with  p_(i,j) = sin(4*PI*i*dx_) + sin(4*PI*j*dy_);
    //grid.valueInitialization(0);

    ///std::cout<<"********* 10 \n";
    ///grid.boundary();

  ///std::cout<<"********* 11 \n";

  while(residual > eps_ && iter <= itermax_)
  {

    ///step 2: SOR solver
    for (int j=1 ; j<=jmax ; j++)
       for (int i=1 ; i<=imax ; i++)
        if (grid.isFluid(i,j)){
        real part1 = coef * ((( grid.p(i+1,j,EAST)+grid.p(i-1,j,WEST))/dx2) + (( grid.p(i,j+1,NORTH)+grid.p(i,j-1,SOUTH))/dy2) - grid.rhs()(i,j));
        grid.p()(i,j) = (1-omg_)*grid.p()(i,j) + part1;
    }

///std::cout<<"********* 12 \n";

    ///step 3: set boundaries
    grid.boundary();

///std::cout<<"********* 13 \n";

    ///step 4: residual calculation
    if (iter % checkfrequency_ == 0){
        residual_calculation(grid);
        residual = resL2Norm(grid);
       // std::cout<<"*************** "<<iter<<"\n";
    }
    iter++;
   // std::cout<<"iteration: "<< iter <<"\tResidule: "<<residual<<"\n";

  }
    if(iter == itermax_) std::cout<<"\nSOR-Calculation is done after ( "<<iter<<" ) iteration with residual norm = "<<residual<<"\n\n";

    if (residual <= eps_)
        return true;
    else return false;
}


void SORSolver::residual_calculation(StaggeredGrid & grid)
{

    int i_max = grid.rhs().getSize(0);
    int j_max =grid.rhs().getSize(1);
    //grid.res().sizing(i_max , j_max);
    real dx2 = grid.dx()* grid.dx();
    real dy2 = grid.dy()*grid.dy();

    for (int j=1 ; j<j_max-1 ; j++)
       for (int i=1 ; i<i_max-1 ; i++)
        if (grid.isFluid(i,j))
          grid.res()(i,j) = grid.rhs()(i,j) - ((( grid.p(i+1,j,EAST)+grid.p(i-1,j,WEST))/dx2) +(( grid.p(i,j+1,NORTH)+grid.p(i,j-1,SOUTH))/dy2)
                            - (2.0)*((dx2+dy2)/(dx2*dy2))*(grid.p()(i,j) ));
}


///L2Norm
real SORSolver::resL2Norm(StaggeredGrid & grid)
{
    int i_max = grid.res().getSize(0);
    int j_max = grid.res().getSize(1);
    real r = 0;
    real residual;

    for (int j=1 ; j<j_max-1 ; j++)
       for (int i=1 ; i<i_max-1 ; i++)
        if (grid.isFluid(i,j))
         r += grid.res()(i,j)*grid.res()(i,j);

    residual = sqrt(r/grid.getNumFluid());

    return residual;
}

// set right hand side values used for testing (with zero value)

void SORSolver::initGridSetup1(StaggeredGrid & grid)
{
    grid.rhs().fill(0);
}

// set right hand side values used for testing (with a sinus function)
void SORSolver::initGridSetup2(StaggeredGrid & grid)
{
    for(int j=1  ; j<grid.rhs().getSize(1)-1; j++)
        for(int i=1  ; i < grid.rhs().getSize(0)-1; i++){
            grid.rhs()(i,j) = sin(2*PI*i*grid.dx());}
}

// Normalize pressure feild in such a way to have a zero avarage
void SORSolver::normalizer(StaggeredGrid & grid)
{
  int imax = grid.imax();
  int jmax = grid.jmax();
  real sum = 0;
  for (int j=1 ; j<=jmax ; j++)
      for (int i=1 ; i<=imax ; i++)
        if (grid.isFluid(i,j))
          sum += grid.p()(i,j);
  real average = sum/grid.getNumFluid();

  for (int j=1 ; j<=jmax ; j++)
     for (int i=1 ; i<=imax ; i++)
        if (grid.isFluid(i,j))
          grid.p()(i,j) = grid.p()(i,j) - average;

}

