#define CHECKFOR_POSNAN
#define CHECKFOR_PROPNAN

//#define SE_CLASS3
//#define CHECKFOR_POSNAN
//#define STOP_ON_ERROR
//#define PRINT_STACKTRACE

double t_stop = 1.0000;

#include "Vector/vector_dist.hpp"
#include "Draw/DrawParticles.hpp"

//#include "doubleToString.h"

int main(int argc, char *argv[])
{
	openfpm_init(&argc, &argv);

	timer timeSim;
	timeSim.start();

	std::string output_dir("/home/ubuntu/zhoulei/github/Pravite-Miscellany/codes/OpenFPM_results/DEM/benchmark/");
	std::string output_filename("DEM_benchmark");

	double m = 1;
	double R = 0.06;
	double cutoff = 2 * R;
	double skin = 0.1 * cutoff;
	constexpr int max_contacts = 36;
	constexpr int max_contacts_def = max_contacts * 3;
	double dt = 0.0001;
	double Im = 2.0 * R * R * m / 5.0;
	double k_n = 7849;
	double k_t = 7849;
	double gamma_n = 3401;
	double gamma_t = 3401;
	double m_eff = 0.5;

	double mu = 0.5;
	double max_vel;
	Vcluster<> &v_cl = create_vcluster();

	Box<3, double> domain({0.0, 0.0, 0.0}, {17.04, 6.0, 6.6});
	Ghost<3, double> g(cutoff + skin);

	//! \cond [constants_def] \endcond

	constexpr int velocity = 0;
	constexpr int force = 1;
	constexpr int omega = 2;
	constexpr int tau = 3;
	constexpr int cpd = 4;
	constexpr int cpi = 5;
	constexpr int ncp = 6;
	constexpr int tt = 7;
	constexpr int gid = 8;

	size_t bc[3] = {NON_PERIODIC, PERIODIC, NON_PERIODIC};

	vector_dist<3, double, aggregate<Point<3, double>, Point<3, double>, Point<3, double>, Point<3, double>, double[max_contacts_def], int[max_contacts], int, int, int>> parts(0, domain, bc, g);

	size_t sz[3] = {143, 51, 56};

	Box<3, double> sand_box({1.8, 0.0, 0.18}, {8.58, 5.9999, 2.7});

	// we draw the initialization
	auto sand_it = DrawParticles::DrawBox(parts, sz, domain, sand_box);

	while (sand_it.isNext())
	{
		// ... add a particle ...
		parts.add();

		// ... and set it position ...
		parts.getLastPos()[0] = sand_it.get().get(0);
		parts.getLastPos()[1] = sand_it.get().get(1);
		parts.getLastPos()[2] = sand_it.get().get(2);

		parts.getLastProp<velocity>()[0] = 0.0;
		parts.getLastProp<velocity>()[1] = 0.0;
		parts.getLastProp<velocity>()[2] = 0.0;

		parts.getLastProp<omega>()[0] = 0.0;
		parts.getLastProp<omega>()[1] = 0.0;
		parts.getLastProp<omega>()[2] = 0.0;

		parts.getLastProp<tau>()[0] = 0.0;
		parts.getLastProp<tau>()[1] = 0.0;
		parts.getLastProp<tau>()[2] = 0.0;

		parts.getLastProp<force>()[0] = 0.0;
		parts.getLastProp<force>()[1] = 0.0;
		parts.getLastProp<force>()[2] = 0.0;

		parts.getLastProp<ncp>() = 0;

		parts.getLastProp<tt>() = 0;

		++sand_it;
	}

	Box<3, double> base_box({0.06, 0.0, 0.06}, {16.98, 5.9999, 0.18});

	// we draw the initialization
	auto base_it = DrawParticles::DrawBox(parts, sz, domain, base_box);

	while (base_it.isNext())
	{
		// ... add a particle ...
		parts.add();

		// ... and set it position ...
		parts.getLastPos()[0] = base_it.get().get(0);
		parts.getLastPos()[1] = base_it.get().get(1);
		parts.getLastPos()[2] = base_it.get().get(2);

		parts.getLastProp<tt>() = 1;

		++base_it;
	}

	Box<3, double> wall_front({16.86, 0.0, 0.06}, {16.98, 5.9999, 6.54});

	// we draw the initialization
	auto wall_f_it = DrawParticles::DrawBox(parts, sz, domain, wall_front);

	while (wall_f_it.isNext())
	{
		// ... add a particle ...
		parts.add();

		// ... and set it position ...
		parts.getLastPos()[0] = wall_f_it.get().get(0);
		parts.getLastPos()[1] = wall_f_it.get().get(1);
		parts.getLastPos()[2] = wall_f_it.get().get(2);

		parts.getLastProp<tt>() = 1;

		++wall_f_it;
	}

	Box<3, double> wall_back({0.06, 0.0, 0.06}, {0.18, 5.9999, 6.54});

	// we draw the initialization
	auto wall_b_it = DrawParticles::DrawBox(parts, sz, domain, wall_back);

	while (wall_b_it.isNext())
	{
		// ... add a particle ...
		parts.add();

		// ... and set it position ...
		parts.getLastPos()[0] = wall_b_it.get().get(0);
		parts.getLastPos()[1] = wall_b_it.get().get(1);
		parts.getLastPos()[2] = wall_b_it.get().get(2);

		parts.getLastProp<tt>() = 1;

		++wall_b_it;
	}

	//! \cond [init sand part] \endcond

	parts.map();
	parts.addComputationCosts(ModelSquare());
	parts.getDecomposition().decompose();
	parts.map();

	// Fill the gid

	auto it_p = parts.getDomainIterator();
	size_t accu = parts.accum();

	while (it_p.isNext())
	{
		auto p = it_p.get();

		parts.getProp<gid>(p) = accu;

		++accu;
		++it_p;
	}

	parts.ghost_get<>();

	size_t cnt = 0;
	size_t cnt_reb = 0;
	auto nlist = parts.getVerlet<VERLETLIST_BAL(3, double)>(cutoff + skin);

	double tot_sp = 0.0;

	double t = 0.0;
	while (t < t_stop)
	{
		auto pit = parts.getDomainIterator();

		max_vel = 0.0;

		//! \cond [integration pos and angvel] \endcond

		// Update
		while (pit.isNext())
		{
			auto p = pit.get();

			if (parts.getProp<tt>(p) == 1)
			{
				++pit;
				continue;
			}

			parts.getProp<velocity>(p)[0] = parts.getProp<velocity>(p)[0] + parts.getProp<force>(p)[0] * dt;
			parts.getProp<velocity>(p)[1] = parts.getProp<velocity>(p)[1] + parts.getProp<force>(p)[1] * dt;
			parts.getProp<velocity>(p)[2] = parts.getProp<velocity>(p)[2] + parts.getProp<force>(p)[2] * dt;

			double norm2 = parts.getProp<velocity>(p)[0] * parts.getProp<velocity>(p)[0] +
										 parts.getProp<velocity>(p)[1] * parts.getProp<velocity>(p)[1] +
										 parts.getProp<velocity>(p)[2] * parts.getProp<velocity>(p)[2];
			if (max_vel < norm2)
			{
				max_vel = norm2;
			}

			parts.getPos(p)[0] = parts.getPos(p)[0] + parts.getProp<velocity>(p)[0] * dt;
			parts.getPos(p)[1] = parts.getPos(p)[1] + parts.getProp<velocity>(p)[1] * dt;
			parts.getPos(p)[2] = parts.getPos(p)[2] + parts.getProp<velocity>(p)[2] * dt;

			if (parts.getPos(p)[0] < domain.getLow(0) || parts.getPos(p)[0] > domain.getHigh(0) ||
					parts.getPos(p)[2] < domain.getLow(2) || parts.getPos(p)[2] > domain.getHigh(2))
			{
				parts.getProp<tt>(p) = 1;
			}

			parts.getProp<omega>(p)[0] = parts.getProp<omega>(p)[0] + parts.getProp<tau>(p)[0] / Im * dt;
			parts.getProp<omega>(p)[1] = parts.getProp<omega>(p)[1] + parts.getProp<tau>(p)[1] / Im * dt;
			parts.getProp<omega>(p)[2] = parts.getProp<omega>(p)[2] + parts.getProp<tau>(p)[2] / Im * dt;

			++pit;
		}
		tot_sp += sqrt(max_vel) * dt;
		v_cl.max(tot_sp);
		v_cl.execute();

		//! \cond [integration pos and angvel] \endcond

		//! \cond [dynamic load balancing] \endcond

		if (tot_sp >= skin / 2.0)
		{
			parts.map();

			// Check if it is time to rebalance

			if (cnt_reb >= 200)
			{
				if (v_cl.rank() == 0)
				{
					std::cout << "Redecomposing" << std::endl;
				}
				cnt_reb = 0;
				parts.addComputationCosts(ModelSquare());
				parts.getDecomposition().redecompose(200);
				parts.map();
			}

			if (v_cl.rank() == 0)
			{
				std::cout << "Reconstruct Verlet" << std::endl;
			}

			parts.ghost_get<velocity, omega, gid>();
			parts.updateVerlet(nlist, cutoff + skin);

			tot_sp = 0.0;
		}
		else
		{
			parts.ghost_get<velocity, omega, gid>();
		}

		//! \cond [dynamic load balancing] \endcond

		//! \cond [calculate force] \endcond

		auto pit2 = parts.getDomainIterator();

		while (pit2.isNext())
		{
			auto p = pit2.get();
			Point<3, double> xp = parts.getPos(p);
			Point<3, double> v_p = parts.getProp<velocity>(p);
			Point<3, double> omega_p = parts.getProp<omega>(p);

			if (parts.getProp<tt>(p) == 1)
			{
				++pit2;
				continue;
			}

			Point<3, double> dF_n({0.0, 0.0, 0.0});
			Point<3, double> dF_t({0.0, 0.0, 0.0});
			Point<3, double> dTau({0.0, 0.0, 0.0});

			int contact_ok[max_contacts];
			for (size_t i = 0; i < max_contacts; i++)
			{
				contact_ok[i] = 0;
			}

			auto NN = nlist.getNNIterator(p.getKey());

			while (NN.isNext())
			{
				auto q = NN.get();

				if (q == p.getKey())
				{
					++NN;
					continue;
				}

				Point<3, double> xq = parts.getPos(q);
				Point<3, double> v_q = parts.getProp<velocity>(q);
				Point<3, double> omega_q = parts.getProp<omega>(q);

				Point<3, double> r_pq = xp - xq;
				double r_s_pq2 = norm2(r_pq);

				// Norm is not defined, next particle
				if (r_s_pq2 == 0)
				{
					continue;
				}

				double delta_ij = 2.0 * R - sqrt(r_s_pq2);

				if (delta_ij < 0.0)
				{
					++NN;
					continue;
				}

				size_t cnt_end = parts.getProp<ncp>(p);
				int this_contact = cnt_end;

				for (size_t k = 0; k < cnt_end; k++)
				{
					if (parts.getProp<cpi>(p)[k] == parts.getProp<gid>(q))
					{
						this_contact = k;
					}
				}

				int cidx;
				if (this_contact == cnt_end)
				{
					parts.getProp<ncp>(p) += 1;
					this_contact = parts.getProp<ncp>(p) - 1;

					cidx = 3 * this_contact;

					if (this_contact >= max_contacts)
					{
						std::cout << "Error reached maximum nunber of contacts points" << std::endl;
					}

					parts.getProp<cpi>(p)[this_contact] = parts.getProp<gid>(q);

					parts.getProp<cpd>(p)[cidx] = 0.0;
					parts.getProp<cpd>(p)[cidx + 1] = 0.0;
					parts.getProp<cpd>(p)[cidx + 2] = 0.0;
				}
				else
				{
					cidx = 3 * this_contact;
				}

				Point<3, double> n_ij = r_pq / sqrt(r_s_pq2);
				Point<3, double> v_rel = v_p - v_q;
				Point<3, double> v_nij = (v_rel * n_ij) * n_ij;
				Point<3, double> v_omega = (omega_p + omega_q) * R;
				Point<3, double> v_cross;

				v_cross.get(0) = v_omega.get(1) * n_ij.get(2) - v_omega.get(2) * n_ij.get(1);
				v_cross.get(1) = v_omega.get(2) * n_ij.get(0) - v_omega.get(0) * n_ij.get(2);
				v_cross.get(2) = v_omega.get(0) * n_ij.get(1) - v_omega.get(1) * n_ij.get(0);

				Point<3, double> v_tij = v_rel - v_nij - v_cross;
				Point<3, double> v_dtij = dt * v_tij;

				parts.getProp<cpd>(p)[cidx] += v_dtij.get(0);
				parts.getProp<cpd>(p)[cidx + 1] += v_dtij.get(1);
				parts.getProp<cpd>(p)[cidx + 2] += v_dtij.get(2);

				Point<3, double> u_ij;

				u_ij.get(0) = parts.getProp<cpd>(p)[cidx];
				u_ij.get(1) = parts.getProp<cpd>(p)[cidx + 1];
				u_ij.get(2) = parts.getProp<cpd>(p)[cidx + 2];

				Point<3, double> F_nij = sqrt(delta_ij / 2 / R) * (k_n * delta_ij * n_ij - gamma_t * m_eff * v_nij);
				dF_n = dF_n + F_nij;

				Point<3, double> F_tij = sqrt(delta_ij / 2 / R) * (-k_t * u_ij - gamma_t * m_eff * v_tij);
				double F_tij_sq = norm2(F_tij);
				double F_nij_sq = mu * mu * norm2(F_nij);
				if (F_tij_sq > F_nij_sq)
				{
					F_tij = F_tij * (F_nij_sq / F_tij_sq);

					parts.getProp<cpd>(p)[cidx] = -1.0 / k_t * (F_tij.get(0) * sqrt(2 * R / delta_ij) + gamma_t * m_eff * v_tij.get(0));
					parts.getProp<cpd>(p)[cidx + 1] = -1.0 / k_t * (F_tij.get(1) * sqrt(2 * R / delta_ij) + gamma_t * m_eff * v_tij.get(1));
					parts.getProp<cpd>(p)[cidx + 2] = -1.0 / k_t * (F_tij.get(2) * sqrt(2 * R / delta_ij) + gamma_t * m_eff * v_tij.get(2));
				}

				dF_t = dF_t + F_tij;
				dTau.get(0) = dTau.get(0) - R * (n_ij.get(1) * F_tij.get(2) - n_ij.get(2) * F_tij.get(1));
				dTau.get(1) = dTau.get(1) - R * (n_ij.get(2) * F_tij.get(0) - n_ij.get(0) * F_tij.get(2));
				dTau.get(2) = dTau.get(2) - R * (n_ij.get(0) * F_tij.get(1) - n_ij.get(1) * F_tij.get(0));

				contact_ok[this_contact] = 1;

				++NN;
			}

			int cnt_end = parts.getProp<ncp>(p);
			int i = 0;
			for (int iread = 0; iread < cnt_end; iread++)
			{
				if (contact_ok[iread] == 1)
				{
					i = i + 1;
					int j = 3 * (i - 1);
					int k = 3 * iread;

					parts.getProp<cpd>(p)[j] = parts.getProp<cpd>(p)[k];
					parts.getProp<cpd>(p)[j + 1] = parts.getProp<cpd>(p)[k + 1];
					parts.getProp<cpd>(p)[j + 2] = parts.getProp<cpd>(p)[k + 2];
				}
			}

			parts.getProp<ncp>(p) = i;

			if (parts.getProp<tt>(p) == 0)
			{
				parts.getProp<force>(p).get(0) = m * 4.905 + dF_n.get(0) + dF_t.get(0);
				parts.getProp<force>(p).get(1) = 0.0 + dF_n.get(1) + dF_t.get(1);
				parts.getProp<force>(p).get(2) = m * -8.49570921 + dF_n.get(2) + dF_t.get(2);

				parts.getProp<tau>(p) = dTau;
			}

			if (parts.getProp<tt>(p) == 1)
			{
				parts.getProp<force>(p) = 0;
				parts.getProp<tau>(p) = 0;
			}

			++pit2;
		}

		//! \cond [calculate force] \endcond

		if (v_cl.rank() == 0)
		{
			std::cout << "simulation TIME: " << t << std::endl;
			std::cout << "clock TIME: " << timeSim.getwct() << std::endl;
		}

		if (cnt % 300 == 0)
		{
			std::cout << "Write " << cnt << std::endl;
			//parts.write_frame("output", cnt);


			parts.write_frame("output",cnt, std::string("time=")+std::to_string(timeSim.getwct()), VTK_WRITER | FORMAT_BINARY);

			//parts.write_frame("output",cnt, VTK_WRITER | FORMAT_BINARY);

			//parts.write_frame("output",cnt, std::string("time=")+std::to_string(timeSim.getwct()));

			
		}

		cnt_reb++;
		cnt++;
		t += dt;
	}

	timeSim.stop();
	if (v_cl.rank() == 0)
	{
		std::cout << "clock TIME: " << timeSim.getwct() << std::endl;
	}

	openfpm_finalize();
}
