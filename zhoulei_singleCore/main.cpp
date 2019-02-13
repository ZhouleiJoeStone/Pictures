#include <iostream>
#include "Vector/vector_dist.hpp"
#include "Draw/DrawParticles.hpp"

#include "SPHKernels.hpp"

#include "SurfaceTension_He2014.hpp"

using namespace SPH;

#include "commom.hpp"

/***********************simulation parameters************************************/
/***********************simulation parameters************************************/

/***********************SPH-related parameters***********************************/

const Real dx = 1e-4; //initial spacing between particles dx 
//const size_t sz[3] = {20,20,20}; // particle number in each direction
const size_t sz[3] = {100,100,100}; // particle number in each direction
const Real dw = 4*dx; //boundary wall width

//!!!Note: the gravity is always in -Y direction for convience in both 2D and 3D simulation
//bounding box 
const Real Xmin = 0.0;
const Real Ymin = 0.0;
const Real Zmin = 0.0;
const Real Xmax = Xmin + sz[0]*dx;
const Real Ymax = Ymin + sz[1]*dx;
const Real Zmax = Zmin + sz[2]*dx;

const Real H = 1.2*dx; //support of the kernel
const Real H_surf = 1.8*dx;//support of the kernel in surface tension 
const Real H_max = std::max(H,H_surf);// cellList cutoff radius
/***********************SPH-related parameters***********************************/



/******************Physical parameters*******************************************/
const Real sigma = 0.5;


// Gravity acceleration
const double gravity = 9.81;

// Reference densitu 1000Kg/m^3
const double rho_zero = 1000.0;

// c_s in the formulas (constant used to calculate the sound speed)
const double coeff_sound = 20.0;

// gamma in the formulas
const double gamma_ = 7.0;

// Maximum height of the fluid water
// is going to be calculated and filled later on
double h_swl = 0.0;

// Filled in initialization
double max_fluid_height = 0.0;

// Eta in the formulas
const double Eta2 = 0.01 * H*H;

// alpha in the formula
const double visco = 0.1;

// cbar in the formula (calculated later)
double cbar = 0.0;

// Filled later require h_swl, it is b in the formulas
double B = 0.0;

// Constant used to define time integration
const double CFLnumber = 0.2;

// Minimum T
const double DtMin = 0.00001;

// Minimum Rho allowed
const double RhoMin = 700.0;

// Maximum Rho allowed
const double RhoMax = 1300.0;

// Mass of the fluid particles
const double MassFluid = 0.000614125;

// Mass of the boundary particles
const double MassBound = 0.000614125;

/******************Physical parameters*******************************************/

/***********************simulation parameters************************************/
/***********************simulation parameters************************************/




#define BOUNDARY 0  // A constant to indicate boundary particles
#define FLUID 1     // A constant to indicate fluid particles



/******************Properties*******************************************************/
const int Type = 0;         //size_t: type of particle :fluid or boundary
const int Mass = 1;         //double: mass 
const int Rho = 2;          //double: density  
const int Rho_prev = 3;     //double: density  at previous step
const int DRho = 4;         //double: Delta rho 
const int Pressure = 5;     //double: pressure
const int Velocity = 6;     //double[DIM]: velocity
const int Velocity_prev = 7;//double[DIM]: velocity at previous step
const int Acc = 8;          //double[DIM]: total acceleration

const int STAcc = 9;        //double[DIM]: acceleration due to surface tension force
const int Color = 10;       //double: color property to caculate surface normal
const int GradC = 11;       //double[DIM]:the gradient of color field
const int GradC2 = 12;      //double: square of the gradient of color field
const int Temperature = 13; //double: temperature 
const int Temperature_pre = 14; //double: temperature 
const int MAcc = 15;        //double[DIM]: acceleration due to tangential Marangoni surface tension force
const int RAcc = 16;        //double[DIM]: acceleration due to recoil pressure

/******************Properties*******************************************************/

// Type of the vector containing particles
typedef vector_dist<DIM,Real,
aggregate<
size_t,     //Type
Real,       //Mass
Real,       //Rho
Real,       //Rho_pre
Real,       //DRho
Real,       //Pressure
Real[DIM],  //Velocity
Real[DIM],  //Velocity_pre
Real[DIM],  //Acc
Real[DIM],  //STAcc
Real,       //Color
Real[DIM],  //GradC
Real,       //GradC2
Real,       //Temperature
Real,       //Temperature_prev
Real[DIM],  //MAcc
Real[DIM]>  //RAcc
> particles;


struct ModelCustom
{
  template<typename Decomposition, typename vector> inline void addComputation( Decomposition & dec,
                                                                                vector & vd,
                                                                                size_t v,
                                                                                size_t p)
  {
    if (vd.template getProp<Type>(p) == FLUID)
      dec.addComputationCost(v,4);
    else
      dec.addComputationCost(v,3);
  }

  template<typename Decomposition> inline void applyModel(Decomposition & dec, size_t v)
  {
    dec.setSubSubDomainComputationCost(v, dec.getSubSubDomainComputationCost(v) * dec.getSubSubDomainComputationCost(v));
  }

  Real distributionTol()
  {
    return 1.01;
  }
};


void drawFluid(particles& vd)
{
  //define fluid vector
  openfpm::vector<Point<DIM,Real>> fluids;

  Real xc = (Xmax+Xmin)/2.0;
  Real yc = (Ymax+Ymin)/2.0;
  Real zc = (Zmax+Zmin)/2.0;

  Real r = (Xmax-Xmin)/4; // bubble radius

  max_fluid_height = 0.0;
  Real min_fluid_height = 99999999;

  if(DIM==3)
  {
    for(Real x = Xmin+dw; x <= Xmax-dw; x +=dx )
      for(Real y = Ymin+dw; y <= Ymax-dw; y +=dx )
        for(Real z = Zmin+dw; z <= Zmax-dw; z +=dx)
        {

          Real xx = x - xc;
          Real yy = y - yc;
          Real zz = z - zc;
          if((xx*xx+yy*yy+zz*zz) > r*r)
          {
            fluids.add({x,y,z});
            max_fluid_height = std::max(max_fluid_height,y);
            min_fluid_height = std::min(min_fluid_height,y);
          } 
          ;
        }
  }

  if(DIM==2)
  {
    for(Real x = Xmin+dw; x <= Xmax-dw; x +=dx )
      for(Real y = Ymin+dw; y <= Ymax-dw; y +=dx )
      {

          Real xx = x - xc;
          Real yy = y - yc;
          if((xx*xx+yy*yy) > r*r)
          {
            fluids.add({x,y});
            max_fluid_height = std::max(max_fluid_height,y);
            min_fluid_height = std::min(min_fluid_height,y);
          } 
      }
  }

  Vcluster & v_cl = create_vcluster();
  v_cl.max(max_fluid_height);
  v_cl.min(min_fluid_height);
  v_cl.execute();


  // here we fill some of the constants needed by the simulation
  h_swl = max_fluid_height - min_fluid_height;
  B = (coeff_sound)*(coeff_sound)*gravity*h_swl*rho_zero / gamma_;
  cbar = coeff_sound * sqrt(gravity * h_swl);


  if (v_cl.getProcessUnitID() == 0)
  {
    std::cout<< "number of fluids: "<< fluids.size() << std::endl;
  }

  // for each particle inside the fluid box ...
  for(size_t i=0; i<fluids.size();i++)
  {
    
    if(vd.getDecomposition().isLocal(fluids.get(i)) == true)
    {
        // Get the position of the probe i
      Point<DIM,Real> rp = fluids.get(i);

      // ... add a particle ...
      vd.add();

      // ... and set it position ...
      vd.getLastPos()[0] = rp.get(0);
      vd.getLastPos()[1] = rp.get(1);
      if(DIM==3)vd.getLastPos()[2] = rp.get(2);

      // and its type.
      vd.template getLastProp<Type>() = FLUID;

      // We also initialize the density of the particle and the hydro-static pressure given by
      //
      // rho_zero*g*h = P
      //
      // rho_p = (P/B + 1)^(1/Gamma) * rho_zero
      //

      vd.template getLastProp<Pressure>() = rho_zero * gravity *  (max_fluid_height - rp.get(1));

      vd.template getLastProp<Rho>() = pow(vd.template getLastProp<Pressure>() / B + 1, 1.0/gamma_) * rho_zero;
      vd.template getLastProp<Rho_prev>() = vd.template getLastProp<Rho>();
      vd.template getLastProp<Velocity>()[0] = 0.0;
      vd.template getLastProp<Velocity>()[1] = 0.0;
      if(DIM==3)vd.template getLastProp<Velocity>()[2] = 0.0;

      vd.template getLastProp<Velocity_prev>()[0] = 0.0;
      vd.template getLastProp<Velocity_prev>()[1] = 0.0;
      if(DIM==3)vd.template getLastProp<Velocity_prev>()[2] = 0.0;

      vd.template getLastProp<Mass>() = MassFluid;
    }
  }
}

void drawWall(particles& vd)
{
  //define fluid vector
  openfpm::vector<Point<DIM,Real>> boundarys;

  if(DIM==3)
  {
    for(Real x = Xmin;x <= Xmax;x +=dx)
      for(Real y = Ymin;y <= Ymax;y +=dx)
        for(Real z = Zmin;z <= Zmax;z +=dx)
        {
          if(x<Xmin+dw||x>Xmax-dw||y<Ymin+dw||y>Ymax-dw||z<Zmin+dw||z>Zmax-dw)
          {
            boundarys.add({x,y,z});
          } 
        }
  }
  else if(DIM==2)
  {
    for(Real x = Xmin;x <= Xmax;x +=dx)
      for(Real y = Ymin;y <= Ymax;y +=dx)
      {
        if(x<Xmin+dw||x>Xmax-dw||y<Ymin+dw||y>Ymax-dw)
        {
          boundarys.add({x,y});
        } 
      }
  }

  Vcluster & v_cl = create_vcluster();
  if (v_cl.getProcessUnitID() == 0)
  {
    std::cout<< "number of boundarys: "<< boundarys.size() << std::endl;
  }

  // for each particle inside the boundary box ...
  for(size_t i=0; i<boundarys.size();i++)
  {
    
    if(vd.getDecomposition().isLocal(boundarys.get(i)) == true)
    {
        // Get the position of the probe i
      Point<DIM,Real> rp = boundarys.get(i);

      // ... add a particle ...
      vd.add();

      // ... and set it position ...
      vd.getLastPos()[0] = rp.get(0);
      vd.getLastPos()[1] = rp.get(1);
      if(DIM==3)vd.getLastPos()[2] = rp.get(2);

      // and its type.
      vd.template getLastProp<Type>() = BOUNDARY;

      vd.template getLastProp<Rho>() = rho_zero;
      vd.template getLastProp<Rho_prev>() = rho_zero;
      vd.template getLastProp<Velocity>()[0] = 0.0;
      vd.template getLastProp<Velocity>()[1] = 0.0;
      if(DIM==3)vd.template getLastProp<Velocity>()[2] = 0.0;

      vd.template getLastProp<Velocity_prev>()[0] = 0.0;
      vd.template getLastProp<Velocity_prev>()[1] = 0.0;
      if(DIM==3)vd.template getLastProp<Velocity_prev>()[2] = 0.0;

      vd.template getLastProp<Mass>() = MassBound;
    }
  }
}


int main(int argc, char ** argv)
{
  
  openfpm_init(&argc,&argv);

  Vcluster & v_cl = create_vcluster();
  timer it_time;

  if (v_cl.getProcessUnitID() == 0)
  std::cout << DIM <<"D simulation " << std::endl;



  #ifdef Dim2
    Box<DIM,Real> domain({Xmin,Ymin},{Xmax,Ymax});//Domain
    size_t bc[DIM] = {NON_PERIODIC,NON_PERIODIC}; //boundary
  #endif
  
  #ifdef Dim3
    Box<DIM,Real> domain({Xmin,Ymin,Zmin},{Xmax,Ymax,Zmax});//Domain
    size_t bc[DIM] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC}; //boundary
  #endif
  
  
  Ghost<DIM,Real> g(4*dx); //ghost layer
  
  particles vd(0,domain,bc,g,DEC_GRAN(512));



  drawFluid(vd);
  drawWall(vd);
  vd.map();

  

  //load balancing
  ModelCustom md;
  vd.addComputationCosts(md);
  vd.getDecomposition().decompose();

  vd.map();

  vd.ghost_get<Type,Mass,Rho,Pressure,Velocity,Color,GradC2,GradC>();

  auto NN = vd.getCellList(2*dx);

  vd.updateCellList(NN);


  WendlandQuinticC2Kernel kernel(H);
  //WendlandQuinticC2Kernel kernel_surface(H_surf);
  WendlandQuinticC2Kernel kernel_surface(H_surf);

  surfaceTensionStep(kernel_surface,vd,NN,sigma);

  vd.map();
  vd.ghost_get<Type,Mass,Rho,Pressure,Velocity,Color,GradC2,GradC>();

  vd.write("Geometry");
    
  

  if (v_cl.getProcessUnitID() == 0)
  {
    std::cout<< "Finished: "<< it_time.getwct() << std::endl;
  }
  
  openfpm_finalize(); 

  return 0;
}

