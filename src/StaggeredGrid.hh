#ifndef STAGGERED_GRID_HH
#define STAGGERED_GRID_HH

#include <iostream>
//#include "Types.hh"
#include "Array.hh"
#include "FileReader.hh"

#include "GrayScaleImage.hh"
#include "lodepng.h"


//*******************************************************************************************************************
/*! Class for storing all arrays required for the simulation
*
* For now it only contains an array for the pressure and another array for the
* right hand side of the pressure equation.
* In following assignments this will be extended and also contain
* arrays for x/y velocity and some arrays for intermediate values.
*
* Feel free to add member functions or variables to this class, but please don't
* delete or rename any existing functions, since future skeletons rely on these functions.
*
*/
//*******************************************************************************************************************
class StaggeredGrid
{
public:
   StaggeredGrid ( ) {}
   // Constructors to manually create staggered grid
   StaggeredGrid ( int xSize, int ySize, real dx, real dy );

   // Constructor to create a staggered grid from a parsed configuration file
   StaggeredGrid ( const FileReader & configuration );

   // Getters / Setters for member variables
   Array & p()    { return p_;    }
   Array & rhs()  { return rhs_;  }
   Array & res()  { return res_;  }
   Array & u()  { return u_;  }
   Array & v()  { return v_;  }
   Array & f()  { return f_;  }
   Array & g()  { return g_;  }
   Array & K()  { return K_;  }
   Array & E()  { return E_;  }
   Array & nueT()  { return nueT_;  } 
   Array & f_mu()  { return f_mu_;  }
   Array & f1()  { return f1_;  }
   Array & f2()  { return f2_;  }
   Array & d()  { return d_;  }
   GrayScaleImage & gray() { return gray_;} // ##

   Array & flagField() { return flagField_; }




   const Array & p()   const { return p_;   }
   const Array & rhs() const { return rhs_; }
   const Array & res() const { return res_; }
   const Array & u() const { return u_;  }
   const Array & v() const { return v_;  }
   const Array & f() const { return f_;  }
   const Array & g() const { return g_;  }
   const Array & K() const { return K_;  }
   const Array & E() const { return E_;  }
   const Array & nueT() const { return nueT_;  }
   const Array & f_mu() const { return f_mu_;  }
   const Array & f1() const { return f1_;  }
   const Array & f2() const { return f2_;  }
   const Array & d() const { return d_;  }
   const GrayScaleImage & gray() const { return gray_;}  // ##


   const Array & flagField() const { return flagField_; }  //if obstacle, it would be 1


   real dx() const { return dx_; }
   real dy() const { return dy_; }
   int imax() const {return imax_; }
   int jmax() const { return jmax_; }

   int xSize() const {return imax_; }
   int ySize() const { return jmax_; }

   StaggeredGrid & operator = (const StaggeredGrid & grid);

   ///initiallize

   void valueInitialization(Array & Arr, real initValue);

   void randomInitialization(Array & Arr);

   void randomInitialization2(Array & Arr);    // initilize with  p_(i,j) = sin(4*PI*i*dx_) + sin(4*PI*j*dy_);
   ///bouandary

   void boundary();

   inline bool isFluid(const int x, const int y);
   inline int getNumFluid();

   // wrapped access
   inline real u(const int x,const int y, Direction dir);
   inline real v(const int x,const int y, Direction dir);
   inline real p(const int x,const int y, Direction dir);
   inline real K(const int x,const int y, Direction dir);
   inline real E(const int x,const int y, Direction dir);

   void createRectangle(real x1, real y1, real x2, real y2);
   void createCircle (real x, real y, real r);
   void setCellToObstacle(int x, int y);

   //PNG functionality
   void writeandSavePNGimage();
   void readPNGimage(const std::string png);

protected:
   Array p_;   //< pressure field
   Array rhs_; //< right hand side of the pressure equation
   Array res_;
   Array v_;
   Array u_;
   Array f_;
   Array g_;
   Array K_;  //Array for turbulent kinetic energy
   Array E_;  //Array for dissipation rate
   Array nueT_;
   Array f_mu_;
   Array f1_;
   Array f2_;
   Array d_;   //contains each field cell distance to obstacles
   Array flagField_;


   real dx_;   //< distance between two grid points in x direction
   real dy_;   //< distance between two grid points in y direction
   int imax_;
   int jmax_;
   GrayScaleImage gray_;

};


inline bool StaggeredGrid::isFluid(const int x, const int y)
{
    if (flagField_(x , y) == 0) return true;
    else return false;
}

inline int StaggeredGrid::getNumFluid()
{
    int sum = 0;
    for(int i=1 ; i<=imax_ ; i++)
        for(int j=1 ; j<=jmax_ ; j++)
          if (isFluid(i,j))
            sum+=1;

    return sum;
}


inline real StaggeredGrid::u(const int x,const int y, Direction dir)
{
    if (isFluid(x,y)==true)
        return u_(x,y);
    else
        switch (dir){
            case NORTH:
                return -u_(x,y-1);
            break;
            case EAST:
                return 0;
            break;
            case SOUTH:
                return -u_(x,y+1);
            break;
            case WEST:
                return 0;
            break;
            default:
                std::cout<< "The wrapped direction of U is wrong! \n";
            return 0;
        }
}

inline real StaggeredGrid::v(const int x,const int y, Direction dir)
{
    if (isFluid(x,y)==true)
        return v_(x,y);
    else
        switch (dir){
            case NORTH:
                return 0;
            break;
            case EAST:
                return -v_(x-1,y);
            break;
            case SOUTH:
                return 0;
            break;
            case WEST:
                return -v_(x+1,y);
            break;
            default:
                std::cout<< "The wrapped direction of V is wrong! \n";
            return 0;
        }
}
inline real StaggeredGrid::p(const int x,const int y, Direction dir)
{
    if (isFluid(x,y)==true)
        return p_(x,y);
    else
        switch (dir){
            case NORTH:
                return p_(x,y-1);
            break;
            case EAST:
                return p_(x-1,y);
            break;
            case SOUTH:
                return p_(x,y+1);
            break;
            case WEST:
                return p_(x+1,y);
            break;
            default:
                std::cout<< "The wrapped direction of P is wrong! \n";
            return 0;
        }

}

inline real StaggeredGrid::K(const int x,const int y, Direction dir)
{
    if (isFluid(x,y)==true)
        return K_(x,y);
    else
        switch (dir){
            case NORTH:
                return -K_(x,y-1);
            break;
            case EAST:
                return -K_(x-1,y);
            break;
            case SOUTH:
                return -K_(x,y+1);
            break;
            case WEST:
                return -K_(x+1,y);
            break;
            default:
                std::cout<< "The wrapped direction of K is wrong! \n";
            return 0;
        }
}

inline real StaggeredGrid::E(const int x,const int y, Direction dir)
{
    if (isFluid(x,y)==true)
        return E_(x,y);
    else
        switch (dir){
            case NORTH:
                return E_(x,y-1);
            break;
            case EAST:
                return E_(x-1,y);
            break;
            case SOUTH:
                return E_(x,y+1);
            break;
            case WEST:
                return E_(x+1,y);
            break;
            default:
                std::cout<< "The wrapped direction of E is wrong! \n";
            return 0;
        }
}


inline void StaggeredGrid::writeandSavePNGimage()
{
  for (int j=jmax_; j!=0; --j)
    for (int i=1; i<=imax_; ++i)
        if(isFluid(i,j))
          gray_.getElement(i-1,jmax_-j)=255;
        else
          gray_.getElement(i-1,jmax_-j)=0;

  for (int j=0; j!=jmax_; ++j)
    if(isFluid(imax_,j))
          gray_.getElement(imax_-1,j)=255;

gray_.save( "Image.png" );

}


//*************************************************************************************************************************
inline void StaggeredGrid::readPNGimage(const std::string pngfile)
{
  GrayScaleImage Gray_(pngfile);
  for (int j=jmax_; j!=0; --j)
    for (int i=1; i!=imax_; ++i)
      if(Gray_.getElement(i-1,jmax_-j)==255)
	flagField_(i,j)=0;
	else flagField_(i,j)=1;

	 Gray_.save( "geometry.png" );
}
#endif //STAGGERED_GRID_HH
