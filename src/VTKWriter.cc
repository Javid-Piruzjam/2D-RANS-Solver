#include "VTKWriter.hh"
#include "Debug.hh"

#include <fstream>
#include <sstream>
#include <iomanip>


template<typename T> struct RealTypeToString         {};
template<>           struct RealTypeToString<float>  { static const char * str; };
template<>           struct RealTypeToString<double> { static const char * str; };

const char * RealTypeToString<float>::str  = "float";
const char * RealTypeToString<double>::str = "double";



VTKWriter::VTKWriter(  const StaggeredGrid & grid, const std::string & basename, bool writePressure, bool writeVelocity, bool writeK, bool writeE )
      : grid_(grid), baseName_( basename ),
        writeVelocity_(writeVelocity), writePressure_(writePressure), writeK_(writeK), writeE_(writeE), counter_ (0 )
{
   ASSERT_MSG( writePressure_ || writeVelocity_ || writeK_ || writeE_ , "VTK Writer has to write at least velocity, pressure, Ka or Epsilon" );

   std::stringstream sstream;
   sstream << "# vtk DataFile Version 4.0\n";
   sstream << "Nusif VTK output\n";
   sstream << "ASCII\n";
   sstream << "DATASET STRUCTURED_POINTS\n";

   sstream << "DIMENSIONS " << grid_.xSize() << " " << grid_.ySize() << " 1\n";
   sstream << "ORIGIN 0 0 0 \n";
   sstream << "SPACING " << grid_.dx() << " " << grid_.dy() << " 1\n";
   sstream << "POINT_DATA " << grid_.xSize() * grid_.ySize() << " \n" << std::endl;

   header_ = sstream.str();
}

void VTKWriter::write()
{
   std::stringstream fileName;
   fileName << baseName_ << "_" <<  std::setw(4) << std::setfill( '0') << counter_ << ".vtk";
   std::ofstream fileStream ( fileName.str().c_str() );

   fileStream << header_;

   if ( writeVelocity_ )
   {
      fileStream << "VECTORS velocity " << RealTypeToString<real>::str << "\n";

      for ( int j = 0; j < grid_.ySize (); ++j )
         for ( int i = 0; i < grid_.xSize (); ++i )
         {
				const real u = 0.5 * ( grid_.u() ( i, j + 1 ) + grid_.u() ( i + 1, j + 1 ) );
				const real v = 0.5 * ( grid_.v() ( i + 1, j ) + grid_.v() ( i + 1, j + 1 ) ) ;

            fileStream << u << " " << v << " " << " 0\n";
         }

      fileStream << "\n";
   }

   if ( writePressure_ )
   {
      fileStream << "SCALARS pressure " << RealTypeToString<real>::str << " 1\n";
      fileStream << "LOOKUP_TABLE default\n";

      for ( int j = 0; j < grid_.ySize (); ++j )
         for ( int i = 0; i < grid_.xSize (); ++i )
            fileStream << grid_.p()( i+1, j+1 ) << "\n";
   }
   
   if ( writeK_ )
   {
      fileStream << "SCALARS K " << RealTypeToString<real>::str << " 1\n";
      fileStream << "LOOKUP_TABLE default\n";

      for ( int j = 0; j < grid_.ySize (); ++j )
         for ( int i = 0; i < grid_.xSize (); ++i )
            fileStream << grid_.K()( i+1, j+1 ) << "\n";
   }
   
   if ( writeE_ )
   {
      fileStream << "SCALARS E " << RealTypeToString<real>::str << " 1\n";
      fileStream << "LOOKUP_TABLE default\n";

      for ( int j = 0; j < grid_.ySize (); ++j )
         for ( int i = 0; i < grid_.xSize (); ++i )
            fileStream << grid_.E()( i+1, j+1 ) << "\n";
   }

   ++counter_;
}
