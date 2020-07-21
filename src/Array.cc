#include "Array.hh"
#include <iostream>
#include <fstream>



//===================================================================================================================
//
//  Constructors
//
//===================================================================================================================


Array::Array( int xSize )
{
   // construct 1D array here
   xSize_= xSize;
   vector_.resize(xSize_,0);
}

Array::Array( int xSize, int ySize )
{
   // construct 2D array here
   xSize_=xSize;
   ySize_=ySize;
   vector_.resize(xSize_*ySize_,0);
}

Array::Array( int xSize, int ySize, int zSize )
{
   // construct 3D array here
   xSize_=xSize;
   ySize_=ySize;
   zSize_=zSize;
   vector_.resize(xSize_*ySize_*zSize_,0);
}

//===================================================================================================================
//
//  Convenience Functions
//
//===================================================================================================================


//initialize the whole array with a constant value
void Array::fill( real value )
{
   /// you might want to use std::fill() here
   std::fill (vector_.begin(), vector_.end(), value);
}


// Print the whole array (for debugging purposes)
void Array::print()
{
   // For 2D Arrays the positive x-coordinate goes to the right
   //                   positive y-coordinate goes upwards
   //      -> the line with highest y-value should be printed first
   std::cout << "\nArray:\n";
   for(int j=ySize_-1; j>=0; j--) {
        for(int i=0; i<xSize_; i++)
            std::cout << vector_[j*xSize_ + i] << "  ";
        std::cout << "\n\n";
   }
   std::cout << '\n';
}


void Array::printToFile()
{
  std::ofstream write;
  write.open("data");
  for(int j=ySize_-2; j>=1; j--){
        for(int i=1; i<xSize_-1; i++)
            write << i << " "<< j << " "<< vector_[j*xSize_ + i]<< "\n";
    write << "\n";
}
}


// return xSize for dimension==0, ySize for dimension==1 and zSize for dimension==2
// other dimension values are not allowed (here it returns a minus value!)
int Array::getSize( int dimension ) const
{
   if (dimension == 0)
        return xSize_;
   else if (dimension == 1)
        return ySize_;
   else if (dimension == 2)
        return zSize_;
   else {
        std::cout<<"wrong dimension! \n";
        return -1;
        }

}

//return total size of the array
int Array::getSize() const
{
   return vector_.size();
}

 void Array::sizing(int xSize, int ySize)
 {
    xSize_ = xSize;
    ySize_ = ySize;
    vector_.resize(xSize*ySize,0);
 }
