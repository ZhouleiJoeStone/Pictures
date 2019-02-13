/*
* SPH kernel functions
*/
#ifndef __SPHKERNELS_hpp__
#define __SPHKERNELS_hpp__

#include "commom.hpp"

#include "Space/Shape/Point.hpp"

namespace SPH
{

	class SPHKernel
	{
	public:
		virtual Real getRadius() = 0;
		virtual void setRadius(Real val) = 0;

		virtual Real W(const Real r) = 0;

		virtual Real W(const Point<DIM,Real> &r)  = 0;

		virtual Point<DIM,Real> gradW(Point<DIM,Real> &r) = 0;
		virtual Real W_zero() = 0;
	};



	/** \brief quintic Wendland C2 kernel.
	$$
	q = r/h
	\\
	W(q) = 
	\begin{cases}
	\alpha_d (1-\frac{q}{2})^4(2q+1), & 0 \le q \le 2 
	\\
	0, & q \ge 2
	\end{cases}

	\\
	\alpha_d = 
	\begin{cases}
	\frac{7}{4 \pi h^{2}}, & DIM =2
	\\
	\frac{21}{16 \pi h^{3}}, & DIM =3
	\end{cases}

	\\

	\nabla_x W(r) = -\frac{1}{h}\alpha_d 5q (1-\frac{q}{2})^3 \frac{dx}{\color{red}{r}} \\
	\nabla_y W(r) = -\frac{1}{h}\alpha_d 5q (1-\frac{q}{2})^3 \frac{dy}{\color{red}{r}} \\
	\nabla_z W(r) = -\frac{1}{h}\alpha_d 5q (1-\frac{q}{2})^3 \frac{dz}{\color{red}{r}} \\

	$$ 
	*/

	class WendlandQuinticC2Kernel
	{
		Real m_radius;
		Real m_k;
		Real m_l;
		Real m_W_zero;
	public:

		WendlandQuinticC2Kernel(Real val) 	{setRadius(val);}

		Real getRadius() { return m_radius; }
		void setRadius(Real val)
		{
			m_radius = val;
			const Real pi = static_cast<Real>(M_PI);

			const Real h3 = m_radius*m_radius*m_radius;
			const Real h5 = h3*m_radius*m_radius;
			if(DIM==3)
			{
				m_k = static_cast<Real>(21.0) / (16.0*pi*h3);
				m_l = -static_cast<Real>(21.0*5.0/ (16.0*pi*h3));
			}
			else if(DIM ==2)
			{
				m_k = static_cast<Real>(7.0) / (4.0*pi*h3);
				m_l = -static_cast<Real>(7.0*5.0/ (4.0*pi*h3));
			}

			m_W_zero = W(0.0);
		}

	public:
		Real W(const Real r)
		{
			Real res = 0.0;
			const Real q = r / m_radius;
			if (q <= 2.0)
			{
				const Real q1= static_cast<Real>(1.0) - q*static_cast<Real>(0.5);
				res = m_k * q1*q1*q1*q1 * (static_cast<Real>(2.0) * q + static_cast<Real>(1.0));
			}
				
			return res;
		}

		Real W(Point<DIM,Real> &r)
		{
			return W(r.norm());
		}

		Point<DIM,Real> gradW(Point<DIM,Real> &r)
		{
			Point<DIM,Real> res;
			const Real rl = r.norm();
			const Real q = rl / m_radius;			
			if (q <= 2.0)
			{
				const Point<DIM,Real> gradq = r / (rl*m_radius);
				const Real q1= static_cast<Real>(1.0) - q*static_cast<Real>(0.5);
				res =  m_l*gradq* q* q1*q1*q1;
			}
			else
				res.zero();

			return res;
		}

		Real W_zero()
		{
			return m_W_zero;
		}
	};


}

#endif