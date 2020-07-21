#include "StaggeredGrid.hh"
#include "Debug.hh"
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265

//class constructor using values
StaggeredGrid::StaggeredGrid( int xSize, int ySize, real dx, real dy ):gray_(xSize, ySize)
{
    imax_= xSize;
    jmax_= ySize;
    CHECK_MSG((imax_>1), "Check the value of imax!\n");
    CHECK_MSG((jmax_>1), "Check the value of jmax!\n");

    p_.sizing(imax_+2 , jmax_+2);
    rhs_.sizing(imax_+2 , jmax_+2);
    res_.sizing(imax_+2 , jmax_+2);
    u_.sizing(imax_+1 , jmax_+2);
    v_.sizing(imax_+2 , jmax_+1);
    f_.sizing(imax_+1 , jmax_+2);
    g_.sizing(imax_+2 , jmax_+1);
    K_.sizing(imax_+2 , jmax_+2);
    E_.sizing(imax_+2 , jmax_+2);
    nueT_.sizing(imax_+2 , jmax_+2);
    f_mu_.sizing(imax_+2 , jmax_+2);
    f1_.sizing(imax_+2 , jmax_+2);
    f2_.sizing(imax_+2 , jmax_+2);
    d_.sizing(imax_+2 , jmax_+2);

    flagField_.sizing(imax_+2 , jmax_+2);
    flagField_.fill(1);
    valueInitialization(flagField(), 0);

    dx_ = dx;
    dy_ = dy;


}

//class constructor using fileReader
StaggeredGrid::StaggeredGrid ( const FileReader & configuration ):gray_(configuration.getIntParameter("imax"),configuration.getIntParameter("jmax"))
{

    //int xSizeT = configuration.getIntParameter("imax") + 2;
    //int ySizeT = configuration.getIntParameter("jmax") + 2;
    imax_= configuration.getIntParameter("imax");
    jmax_= configuration.getIntParameter("jmax");

    CHECK_MSG((imax_>1), "Check the value of imax!\n");
    CHECK_MSG((jmax_>1), "Check the value of jmax!\n");

    p_.sizing(imax_+2 , jmax_+2);
    rhs_.sizing(imax_+2 , jmax_+2);
    res_.sizing(imax_+2 , jmax_+2);
    u_.sizing(imax_+1 , jmax_+2);
    v_.sizing(imax_+2 , jmax_+1);
    f_.sizing(imax_+1 , jmax_+2);
    g_.sizing(imax_+2 , jmax_+1);
    K_.sizing(imax_+2 , jmax_+2);
    E_.sizing(imax_+2 , jmax_+2);
    nueT_.sizing(imax_+2 , jmax_+2);
    f_mu_.sizing(imax_+2 , jmax_+2);
    f1_.sizing(imax_+2 , jmax_+2);
    f2_.sizing(imax_+2 , jmax_+2);
    d_.sizing(imax_+2 , jmax_+2);

    flagField_.sizing(imax_+2 , jmax_+2);
    flagField_.fill(1);
    valueInitialization(flagField(), 0);

    dx_ = configuration.getRealParameter("xlength") / (configuration.getIntParameter("imax"));
    dy_ = configuration.getRealParameter("ylength") / (configuration.getIntParameter("jmax"));


}

 StaggeredGrid & StaggeredGrid::operator = (const StaggeredGrid & grid)
 {
    p_=grid.p();   //< pressure field
    rhs_=grid.rhs(); //< right hand side of the pressure equation
    res_=grid.res();
    v_=grid.v();
    u_=grid.u();
    f_=grid.f();
    g_=grid.g();

    dx_=grid.dx();   //< distance between two grid points in x direction
    dy_=grid.dy();   //< distance between two grid points in y direction
    imax_=grid.imax();
    jmax_=grid.jmax();
    return *this;
 }

///initiallize pressure field with a constant value
void StaggeredGrid::valueInitialization(Array & Arr, real initValue)
{
    int iSize = Arr.getSize(0);
    int jSize = Arr.getSize(1);
    for (int j=1 ; j<jSize-1 ; j++){
        for (int i=0 ; i<=iSize-1 ; i++){
            if (isFluid(i,j))
              Arr(i,j) = initValue;
            else
              Arr(i,j) = 0;
        }
  }

}

///initiallize pressure field with random values
void StaggeredGrid::randomInitialization(Array & Arr)
{
    srand (time(NULL));
    int iSize = Arr.getSize(0);
    int jSize = Arr.getSize(1);
    for (int j=1 ; j<jSize-1 ; j++)
        for (int i=1 ; i<iSize-1 ; i++){
            if (isFluid(i,j))
              Arr(i,j) = ((real) rand() / (RAND_MAX));
            else
              Arr(i,j) = 0;
        }
}

///initiallize pressure field with a function
void StaggeredGrid::randomInitialization2(Array & Arr)   // initilize with  p_(i,j) = sin(4*PI*i*dx_) + sin(4*PI*j*dy_);
{
    srand (time(NULL));
    int iSize = Arr.getSize(0);
    int jSize = Arr.getSize(1);
    for (int j=1 ; j<jSize-1 ; j++)
        for (int i=1 ; i<iSize-1 ; i++){
            if (isFluid(i,j))
              Arr(i,j) = sin(4*PI*i*dx_) + sin(4*PI*j*dy_);
            else
              Arr(i,j) = 0;
        }
}

///bouandary
void StaggeredGrid::boundary()
{
    int imax = imax_;
    int jmax = jmax_;


    for (int i=1 ; i<=imax ; i++) {
        p_(i,0) = p_(i,1);
        p_(i,jmax+1) = p_(i,jmax);

    }
    for (int j=1 ; j<=jmax ; j++){
        p_(0,j) = p_(1,j);
        p_(imax+1, j) = p_(imax, j);
    }

}

void StaggeredGrid::createRectangle(real x1, real y1, real x2, real y2)
{
    int X1 = x1/dx_;
    int Y1 = y1/dy_;
    int X2 = x2/dx_;
    int Y2 = y2/dy_;

    for (int j=Y1+1 ; j<=Y2 ; j++)
        for(int i=X1+1 ; i<=X2 ; i++)
          setCellToObstacle(i,j);
}
void StaggeredGrid::createCircle (real x, real y, real r)
{
    int X1 = (x-r)/dx_;
    int Y1 = (y-r)/dy_;
    int X2 = (x+r)/dx_;
    int Y2 = (y+r)/dy_;

    for (int j=Y1 ; j<=Y2 ; j++)
        for(int i=X1 ; i<=X2 ; i++)
          if ( (i*dx_-x)*(i*dx_-x) + (j*dy_-y)*(j*dy_-y) <= r*r   )
            setCellToObstacle(i,j);
}

void StaggeredGrid::setCellToObstacle(int x, int y)
{
    flagField_(x,y)=1;
}
