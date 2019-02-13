#include <iostream>
#include "Vector/vector_dist.hpp"

#include "SPHKernels.hpp"

const int x = 0; //x坐标
const int y = 1; //y坐标

//constexpr int A = 0; //属性1
//constexpr int B = 1; //属性2

const int A = 0; //属性1
const int B = 1; //属性2


#define H 2.0

const double alpha_d = 21.0/(16.0*M_PI*H*H*H);

inline double Wab(double r)
{
   double q = r/H;
   double q05 = 1-q*0.5;

   if (q < 2)
      return alpha_d*q05*q05*q05*q05*(2*q+1);
   else
      return 0.0;
}


inline void DWab(Point<3,double> & dx, Point<3,double> & DW, double r, bool print)
{
   double factor = 0;

   const double q=r/H;
   if(q<2)
   {
      const double q05 = 1-q*0.5;
      factor = -alpha_d*5*q*q05*q05*q05 / H /r;
   }

    DW.get(0) = factor * dx.get(0);
    DW.get(1) = factor * dx.get(1);
    DW.get(2) = factor * dx.get(2);
}



int main(int argc, char ** argv)
{
    openfpm_init(&argc,&argv);
    
    Box<2,float> domain({0.0, 0.0},{1.0, 1.0});//Domain
    size_t bc[2] = {PERIODIC,PERIODIC}; //boundary
    Ghost<2,float> g(0.001); //ghost layer
    
    vector_dist<2,float,aggregate<float,float> > vd(10,domain,bc,g);
    
    for(size_t i = 0; i<vd.size_local(); i++)
    {
        vd.getPos(i)[x] = (float)rand()/RAND_MAX;
        vd.getPos(i)[y] = (float)rand()/RAND_MAX;
        
        vd.getProp<A>(i) = 0.5;
        vd.getProp<B>(i) = 0.6;
    }
    
    std::cout<< "Hello world"<<std::endl;


    WendlandQuinticC2Kernel<3,double> WC2(2.0);
    std::cout<< "W_zero(): " << WC2.W_zero() << std::endl;
    std::cout<< "Wab(0): " << Wab(0)<< std::endl;


    std::cout<< "W(1)): " << WC2.W(1)<< std::endl;
    std::cout<< "Wab(1): " << Wab(1)<< std::endl;

    Point<3,double> p({0.1,0.1,0.2});

    double r = p.norm();
    Point<3,double> DW;
    DWab(p,DW,r,false);

    std::cout<< "gradW(p): "<< WC2.gradW(p)[0] <<" "<<WC2.gradW(p)[1] <<" "<<WC2.gradW(p)[2] <<" "<<std::endl;
    std::cout<< "DWab(p): "<< DW.get(0) <<" "<<DW[1] <<" "<<DW[2] <<" "<<std::endl;

    


    
    openfpm_finalize();
    
    return 0;
}

