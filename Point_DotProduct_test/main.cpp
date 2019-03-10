
#include <stddef.h>
#include "Vector/vector_dist.hpp"

typedef double Real;

int main(int argc, char* argv[])
{
    // initialize the library
	openfpm_init(&argc,&argv);

	Point<2,Real> p0({0.1,0.1});
	Point<2,Real> p1({0.1,0.1});
	
	Real p01 = p0*p1;
	Real p01_ = (p0*p1);
	std::cout<< "p0*p1   : "<< p01 <<"\n";
	std::cout<< "(p0*p1)   : "<< p01_ <<"\n";

	Point<2,Real> p02 = 2*p0;
	Point<2,Real> p020 = 2.0*p0;
	std::cout<< "2*p0   : "<< p02[0] <<"  " << p02[1] <<"\n";
	std::cout<< "2.0*p0 : "<< p020[0] <<"  " << p020[1] <<"\n";


	Real p2 = 2*(p0*p1);
	Real p3 = 2.0*(p0*p1);
	std::cout<< "2*(p0*p1): "<< p2 <<"\n";
	std::cout<< "2.0*(p0*p1): "<< p3 <<"\n";

	Point<2,Real> pp = (p0*p1)*p1;
	std::cout<< "(p0*p1)*p1: "<< pp[0] <<" "<<pp[1] <<"\n";

	Point<2,Real> p22 = 2*(p0*p1)*p1;
	Point<2,Real> p33 = 2.0*(p0*p1)*p1;
	std::cout<< "2*(p0*p1)*p1: "<< p22[0] <<" "<<p22[1] <<"\n";
	std::cout<< "2.0*(p0*p1)*p1: "<< p33[0]<<" "<<p33[1] <<"\n";


	openfpm_finalize();
}
