#include "FluidSimulator.hh"
#include "Debug.hh"
#include <iostream>
#include <math.h>

int main()
{
  FileReader reader;

  PROG("Registering the Parameters...");

  reader.registerStringParameter ("boundary_condition_N");
  reader.registerStringParameter ("boundary_condition_S");
  reader.registerStringParameter ("boundary_condition_E");
  reader.registerStringParameter ("boundary_condition_W");
  reader.registerRealParameter ("boundary_velocity_N");
  reader.registerRealParameter ("boundary_velocity_S");
  reader.registerRealParameter ("boundary_velocity_E");
  reader.registerRealParameter ("boundary_velocity_W");
  reader.registerRealParameter ("GX");
  reader.registerRealParameter ("GY");
  reader.registerRealParameter ("Re");
  reader.registerRealParameter ("U_INIT");
  reader.registerRealParameter ("V_INIT");
  reader.registerRealParameter ("P_INIT");
  reader.registerRealParameter("xlength");
  reader.registerRealParameter("ylength");
  reader.registerIntParameter("imax");
  reader.registerIntParameter("jmax");
  reader.registerRealParameter ("RectangleX1");
  reader.registerRealParameter ("RectangleY1");
  reader.registerRealParameter ("RectangleX2");
  reader.registerRealParameter ("RectangleY2");
  reader.registerRealParameter ("CircleX");
  reader.registerRealParameter ("CircleY");
  reader.registerRealParameter ("CircleR");
  
  
  reader.registerRealParameter ("K_INIT");
  reader.registerRealParameter ("E_INIT");
  reader.registerRealParameter ("C_e");
  reader.registerRealParameter ("C_mu");
  reader.registerRealParameter ("C1");
  reader.registerRealParameter ("C2");

  reader.registerRealParameter ("dt");
  reader.registerIntParameter ("timesteps");
  reader.registerRealParameter("safetyfactor");
  reader.registerIntParameter ("itermax");
  reader.registerRealParameter ("eps");
  reader.registerRealParameter ("omg");
  reader.registerRealParameter ("gamma");
  reader.registerIntParameter ("checkfrequency");
  reader.registerIntParameter ("normalizationfrequency");
  reader.registerIntParameter ("outputinterval");


  reader.readFile("backstep.par");
  FluidSimulator simul (reader);
  simul.grid().createRectangle(reader.getRealParameter("RectangleX1"),reader.getRealParameter("RectangleY1"),reader.getRealParameter("RectangleX2"),reader.getRealParameter("RectangleY2"));
  //simul.grid().createCircle(reader.getRealParameter("CircleX"),reader.getRealParameter("CircleY"),reader.getRealParameter("CircleR"));
  simul.grid().writeandSavePNGimage();

  //simul.simulate(20);
  //simul.simulateTimeStepCount(reader.getIntParameter("timesteps"));
  simul.simulateTurbulence(reader.getIntParameter("timesteps"));
  

  return 0;

}
