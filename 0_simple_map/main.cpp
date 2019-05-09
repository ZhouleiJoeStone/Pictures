#include <stddef.h>
#include "Vector/vector_dist.hpp"

int main(int argc, char *argv[])
{

	openfpm_init(&argc, &argv);

	Box<2, float> domain({0.0, 0.0}, {10e6, 10e6});
	//Box<2, float> domain({0.0, 0.0}, {1.0, 1.0});

	size_t bc[2] = {PERIODIC, PERIODIC};

	Ghost<2, float> g(0.01);

	vector_dist<2, float, aggregate<float>> vd(0, domain, bc, g);

	Vcluster<> &v_cl = create_vcluster();

	if (v_cl.getProcessUnitID() == 0)
	{
		for (size_t i = 0; i < 1000000; i++)
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

	vd.write_frame("particles", 0);

	for (int i = 1; i < 10; i++)
	{

		//move particles
		auto it = vd.getDomainIterator();
		while (it.isNext())
		{
			auto key = it.get();
			vd.getPos(key)[0] += 1 ;
			vd.getPos(key)[1] += 1 ;
			++it;
		}

		vd.map();
		vd.addComputationCosts();
		vd.getDecomposition().decompose();
		vd.map();

		vd.ghost_get<>();
		vd.write_frame("particles", i);
	}

	openfpm_finalize();
}
