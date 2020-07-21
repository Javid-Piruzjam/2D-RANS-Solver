#ifndef ARRAY_HH
#define ARRAY_HH

#include "Types.hh"
#include <vector>
#include <assert.h>


//*******************************************************************************************************************
/*!  Array class for 1,2 and 3 dimensions
*
*    - all elements should be stored in a contiguous chunk of memory ( no vector<vector> ! )
*/
//*******************************************************************************************************************
class Array
{
public:
   // Constructors for 1D,2D and 3D
   Array() {vector_.resize(1,0);}
   Array( int xSize );
   Array( int xSize, int ySize );
   Array( int xSize, int ySize, int zSize );


   // Depending on your implementation you might need the following:
   /// ~Array(); it seems to be unnecessary, as it is usually needed in the case of using pointer stuff!
   // Array(const Array& s);
   // Array& operator= (const Array& s);

   // Access Operators for 1D, 2D and 3D
   inline real & operator () ( int i );
   inline real & operator () ( int i ,int j );
   inline real & operator () ( int i, int j, int k );

   // for const Arrays the following access operators are required
   inline const real & operator () ( int i ) const;
   inline const real & operator () ( int i ,int j ) const;
   inline const real & operator () ( int i, int j, int k ) const;

   // initialize the whole array with a constant value
   void fill( real value );

   // return total size of the array
   int getSize() const;

   // return xSize for dimension==0, ySize for dimension==1 and zSize for dimension==2
   // other dimension values are not allowed
   int getSize(int dimension ) const;

   ///
   void sizing (int xSize, int ySize);

   // Print the whole array ( for debugging purposes )
   void print();

   void printToFile();

   std::vector <real> contain() const {return vector_;}

private:

    int xSize_;
    int ySize_;
    int zSize_;
    std::vector <real> vector_;

};


//===================================================================================================================
//
//  Inline Access Operators and Sizes
//
//===================================================================================================================


// Operator() 1D
inline real& Array::operator ()(int i)
{
   assert(i<xSize_ && i>=0);
   return vector_[i];
}

// Operator() 2D
inline real& Array::operator ()(int i,int j)
{

   assert(i>=0 && j>=0 && i<xSize_ && j<ySize_);
   return vector_[j*xSize_ + i];
}

// Operator() 3D
inline real& Array::operator ()(int i, int j, int k)
{
   assert(i>=0 && j>=0 && k>=0 && i<xSize_ && j<ySize_ && k<zSize_);
   return vector_[k*xSize_*ySize_ + j*xSize_ + i];
}

inline const real & Array::operator () ( int i ) const
{
   assert(i<xSize_ && i>=0);
   return vector_[i];
}

inline const real & Array::operator () ( int i ,int j ) const
{
   assert(i>=0 && j>=0 && i<xSize_ && j<ySize_);
   return vector_[j*xSize_ + i];
}
inline const real & Array::operator () ( int i, int j, int k ) const
{
   assert(i>=0 && j>=0 && k>=0 && i<xSize_ && j<ySize_ && k<zSize_);
   return vector_[k*xSize_*ySize_ + j*xSize_ + i];
}


#endif //ARRAY_HH

