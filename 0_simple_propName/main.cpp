
#include <stddef.h>
#include "Vector/vector_dist.hpp"

constexpr int gid = 0;				 //global id
constexpr int material_id = 1; //材料索引
constexpr int radius = 2;
constexpr int mass = 3;

constexpr int force = 4;
constexpr int tau = 5;
constexpr int velocity = 6;
constexpr int omega = 7;

constexpr int ncp = 8;			 //number of contact points
constexpr int cpi = 9;			 //particle index of the contact point
constexpr int cpd = 10;			 //contact point deformation
constexpr int radius_e = 11; // equivalent radius
constexpr int mass_e = 12;	 //equivalent mass

//constexpr int max_contacts = 3;
constexpr int max_contacts = 4;

typedef vector_dist<3, double,
										aggregate<
												int,										 //0 gid: global id unique for each particles
												int,										 //1 material_id:
												double,									 //2 radius
												double,									 //3 mass
												Point<3, double>,				 //4 force:	particle force
												Point<3, double>,				 //5 tau:		moment of the force
												Point<3, double>,				 //6 velocity: particle velocity
												Point<3, double>,				 //7 omega:	angular velocity
												int,										 //8 ncp: number of contact points
												int[max_contacts],			 //9 cpi: 		particle index of the contact point
												double[max_contacts][3], //10 cpd: contact point deformation
												double[max_contacts],		 //11 radius_e: equivalent radius
												double[max_contacts]		 //12 mass_e: equivalent mass
												>>
		particles;

int main(int argc, char *argv[])
{

	openfpm_init(&argc, &argv);

	Box<3, float> domain({0.0, 0.0,0.0}, {1.0, 1.0,1.0});
	size_t bc[3] = {PERIODIC, PERIODIC, PERIODIC};

	Ghost<3, float> g(0.01);

	particles vd(1, domain, bc, g);

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		// we define x, assign a random position between 0.0 and 1.0
		vd.getPos(key)[0] = (float)rand() / RAND_MAX;

		// we define y, assign a random position between 0.0 and 1.0
		vd.getPos(key)[1] = (float)rand() / RAND_MAX;

		// next particle
		++it;
	}

	vd.map();

	openfpm::vector<std::string> particlePropNames({"gid",
																									"material_id",
																									"radius",
																									"mass",
																									"force",
																									"tau",
																									"velocity",
																									"omega",
																									"ncp",
																									"cpi",
																									"cpd",
																									"radiusEffective",
																									"massEffective"});
	vd.setPropNames(particlePropNames);

	// save vtk format (vtk is always the default)
	vd.write("particles");

	openfpm_finalize();
}
