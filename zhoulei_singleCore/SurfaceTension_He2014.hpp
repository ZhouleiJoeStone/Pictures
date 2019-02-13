#ifndef __SurfaceTension_He2014_h__
#define __SurfaceTension_He2014_h__


#include "commom.hpp"


#include "Vector/vector_dist.hpp"
#include "SPHKernels.hpp"

#define BOUNDARY 0  // A constant to indicate boundary particles
#define FLUID 1     // A constant to indicate fluid particles

extern const int Type;//size_t: type of particle :fluid or boundary
extern const int Mass; //mass property
extern const int Rho; //density property
extern const int Color; // color property
extern const int GradC;       //double[DIM]:the gradient of color field
extern const int GradC2;      //square of the gradient of color field
extern const int STAcc; // acceleration due to surface tension force

template<typename Kernel, typename Vector, typename CellList>
void surfaceTensionStep(Kernel & kernel,Vector & vd,CellList & NN, Real sigma)
{
	computeColorField(kernel,vd,NN);
	computeGradientColorField(kernel,vd,NN);
	computeSurfaceTensionForce(kernel,vd,NN,sigma);
}

template<typename Kernel, typename Vector, typename CellList>
void computeColorField(Kernel & kernel,Vector & vd,CellList & NN)
{
	// Get the iterator
    auto it = vd.getDomainIterator();
    // For each particle ...

    const Real H2 = 2.0*kernel.getRadius();

    while (it.isNext())
    {
        // ... a
        auto a = it.get();
        // Get the color of the particle a
        vd.template getProp<Color>(a) = 0;

        vd.template getProp<Color>(a) += vd.template getProp<Mass>(a) / vd.template getProp<Rho>(a) * kernel.W_zero();


        // Get the position of the particle a
        Point<DIM,Real> xa = vd.template getPos(a);

        // Get an iterator of all the particles neighborhood of a
        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
        // For each particle near a
        while (Np.isNext())
        {
            // Get the particle b near to a
            auto b = Np.get();
            // if (b == a) skip this particle
            if (a.getKey() == b)    {++Np; continue;};

            // Get the position of the particle b
            Point<DIM,Real> xb = vd.getPos(b);
            // Calculate the distance vector between a and b
            Point<DIM,Real> dr = (xa - xb);
            // Calculate the distance between a and b
            Real r = dr.norm();

            if(r< H2)
            {

        		vd.template getProp<Color>(a) += vd.template getProp<Mass>(b) / vd.template getProp<Rho>(b) * kernel.W(r);
            }
   
            ++Np;
        }
        // Next particle p
        ++it;
    }
}

template<typename Kernel, typename Vector, typename CellList>
void computeGradientColorField(Kernel & kernel,Vector & vd,CellList & NN)
{
	// Get the iterator
    auto it = vd.getDomainIterator();
    // For each particle ...

    const Real H2 = 2.0*kernel.getRadius();

    while (it.isNext())
    {
        // ... a
        auto a = it.get();
        // Get the color gradient of the particle a
        Point<DIM,Real> gradCa;
        gradCa.zero();

        //Don't add the color gradient at particle a ,because gradW is equal to 0
        //gradCa += vd.template getProp<Mass>(a) / vd.template getProp<Rho>(a) *vd.template getProp<Color>(a) * kernel.gradW(Point<DIM,Real>({0,0,0}));

        // Get the position of the particle a
        Point<DIM,Real> xa = vd.template getPos(a);

        // Get an iterator of all the particles neighborhood of a
        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
        // For each particle near a
        while (Np.isNext())
        {
            // Get the particle b near to a
            auto b = Np.get();

            // if (b == a) skip this particle
            if (a.getKey() == b)    {++Np; continue;};

            // Get the position of the particle b
            Point<DIM,Real> xb = vd.getPos(b);
            // Calculate the distance vector between a and b
            Point<DIM,Real> dr = (xa - xb);
            // Calculate the distance between a and b
            Real r = dr.norm();

            if(r < H2)
            {

        		gradCa += vd.template getProp<Mass>(b) / vd.template getProp<Rho>(b) *vd.template getProp<Color>(b) * kernel.gradW(dr);
            }
   
            ++Np;
        }

        gradCa = gradCa* (static_cast<Real>(1.0) / vd.template getProp<Color>(a));

        
        vd.template getProp<GradC>(a)[0] =  gradCa.get(0);
        vd.template getProp<GradC>(a)[1] =  gradCa.get(1);
        if(DIM==3)  vd.template getProp<GradC>(a)[2] =  gradCa.get(2); 

        vd.template getProp<GradC2>(a) = norm2(gradCa);

        // Next particle p
        ++it;
    }
}

template<typename Kernel, typename Vector, typename CellList>
void computeSurfaceTensionForce(Kernel & kernel,Vector & vd,CellList & NN, Real sigma)
{
    // Get the iterator
    auto it = vd.getDomainIterator();
    // For each particle ...

    const Real H2 = 2.0*kernel.getRadius();

    while (it.isNext())
    {
        // ... a
        auto a = it.get();

        Real gradC2a = vd.template getProp<GradC2>(a);
        Real rhoa = vd.template getProp<Rho>(a);
        Real factor  = static_cast<Real>(0.25)*sigma/ rhoa;
        Point<DIM,Real> stAcc;
        stAcc.zero();

        if(vd.template getProp<Type>(a) == FLUID)
        {
            // Get the position of the particle a
            Point<DIM,Real> xa = vd.template getPos(a);

            // Get an iterator of all the particles neighborhood of a
            auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
            // For each particle near a
            while (Np.isNext())
            {
                // Get the particle b near to a
                auto b = Np.get();
                // if (b == a) skip this particle
                if (a.getKey() == b)    {++Np; continue;};

                // Get the position of the particle b
                Point<DIM,Real> xb = vd.getPos(b);
                // Calculate the distance vector between a and b
                Point<DIM,Real> dr = (xa - xb);
                // Calculate the distance between a and b
                Real r = dr.norm();

                if(r < H2)
                {
                    stAcc += factor*vd.template getProp<Mass>(b)/vd.template getProp<Rho>(b)*(vd.template getProp<GradC2>(a) + vd.template getProp<GradC2>(b))* kernel.gradW(dr);
                }
                ++Np;
            }

        }

        

        vd.template getProp<STAcc>(a)[0] =  stAcc.get(0);
        vd.template getProp<STAcc>(a)[1] =  stAcc.get(1);
        if(DIM==3)  vd.template getProp<STAcc>(a)[2] =  stAcc.get(2); 

        // Next particle p
        ++it;
    }
}

#endif