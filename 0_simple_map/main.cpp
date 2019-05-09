#include <stddef.h>
#include "Vector/vector_dist.hpp"

int main(int argc, char *argv[])
{

	openfpm_init(&argc, &argv);

	Box<2, float> domain({0.0, 0.0}, {10.0, 10.0});
	//Box<2, float> domain({0.0, 0.0}, {1.0, 1.0});

	size_t bc[2] = {PERIODIC, PERIODIC};

	Ghost<2, float> g(0.01);

	vector_dist<2, float, aggregate<float>> vd(0, domain, bc, g);

	Vcluster<> &v_cl = create_vcluster();

	if (v_cl.getProcessUnitID() == 0)
	{
		for (size_t i = 0; i < 100; i++)
		{
			vd.add();
			// we define x, assign a random position between 0.0 and 1.0
			vd.getLastPos()[0] = (float)rand() / RAND_MAX;
			// we define y, assign a random position between 0.0 and 1.0
			vd.getLastPos()[1] = (float)rand() / RAND_MAX;
		}
	}

	vd.map();
	vd.addComputationCosts();
	vd.getDecomposition().decompose();
	vd.map();

	vd.ghost_get<>();
	vd.write("particles");

	openfpm_finalize();
}
