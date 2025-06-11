#include "cgle_sim.h"

#include <iostream>
#include <string>
#include <chrono>
#include <random>
#include <numer_visualize.h>
#include <numer_fourier.h>
#include <numer_grid.h>
#include <numer_file.h>

#include <stbi.h>





using namespace numer;

static const char* WF_FILE_PATH = "./image/simul/cgle/glfield";
constexpr unsigned FREAM_INTERVAL = 5;
constexpr double dt = 1.0 / (1.0 * FREAM_INTERVAL);
constexpr unsigned N1 = 720;
constexpr unsigned N2 = 1280;
constexpr double L1 = 500.0;



//∂Ψ/∂t = (ε + ib)∇²Ψ + cΨ - (d + ie)|Ψ|²Ψ
constexpr double ep = 1.0;
constexpr double b = 0.0;
constexpr double c = 1.0;
constexpr double d = 1.0;
constexpr double e = 1.2;

static mat<Complex> initialize(unsigned N1_, unsigned N2_) {
	mat<Complex> result(N1_, N2_);

	unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine rand_eng(seed);
	std::uniform_real_distribution<double> amp_distr(0.8, 1.0);
	std::uniform_real_distribution<double> pha_distr(0.0, 2.0 * Pi);

	result.refill([&](size_t i, size_t j) {
		if (j <= N2_ / 2)
			return Complex::expi(pha_distr(rand_eng)) * amp_distr(rand_eng);
		else
			return Complex::zero();
		});

	return result;
}


CGLESim::CGLESim(unsigned Steps)
	: Sim(Steps), field_{}
{

	{//initialization logic
		unsigned op = 0;
	init:
		std::cout << "\n\tchoose an initialization:\n\t\t0.new computation\n\t\t1.load from file\n\n\tenter: ";
		std::cin >> op;

		switch (op) {
		case 0:
			profile_.serial = 0;
			profile_.frame_interval = FREAM_INTERVAL;
			profile_.dt = dt;
			profile_.N1 = N1;
			profile_.N2 = N2;
			profile_.L1 = L1;
			profile_.ep = ep;
			profile_.b = b;
			profile_.c = c;
			profile_.d = d;
			profile_.e = e;
			field_ = std::move(initialize(N1, N2));
			break;
		case 1:
		tryagain:
			if (readMat_strong(WF_FILE_PATH, field_, profile_)) break;
			else {
				std::cout << "\n\t\"glfield.nmmat\" doesn't exist or isn't usable !\n\t\t0.try again\n\t\t1.go back\n\n\tenter:";
				unsigned op1 = 0;
				std::cin >> op1;
				switch (op1) {
				case 0: goto tryagain;
				case 1: goto init;
				default: goto init;
				}
			}
		default: goto init;
		}
	}//end initialization logic

	std::cout << std::endl;
}





mat<RGB> CGLESim::renderFrame() {

	//----------------------------------time parameter---------------------------------

	//fixed step length

	//------------------------------------preparation----------------------------------

	//∂Ψ/∂t = (ε + ib)∇²Ψ + cΨ - (d + ie)|Ψ|²Ψ

	const double L2 = profile_.L1 * profile_.N2 / profile_.N1;

	const auto nonlin_phs = [=](const Complex& spatio_) {
		constexpr Complex i = Complex::i();
		return exp(-(profile_.d + i * profile_.e) * norm(spatio_) * profile_.dt / 2.0);
		};

	const FFTFreq freq1(profile_.N1, 2.0 * profile_.L1);
	const FFTFreq freq2(profile_.N2, 2.0 * L2);

	const auto linear_phs = [=](size_t i_, size_t j_) {
		constexpr Complex i = Complex::i();
		double kx = freq1(i_), ky = freq2(j_);
		double k2 = kx * kx + ky * ky;
		return exp((profile_.c - (profile_.ep + i * profile_.b) * k2) * profile_.dt);
		};



	//-------------------------------------simulation----------------------------------

	thread_local mat<Complex> nlphs_factors(profile_.N1, profile_.N2);
	thread_local auto lphs_factors = mat<Complex>::creat_par(profile_.N1, profile_.N2, linear_phs);

	const auto tssp1 = [&]() {
		nlphs_factors.refill_par(field_, nonlin_phs);
		field_.overlay_par(nlphs_factors, multiply<Complex, Complex>);
		};

	const auto tssp2 = [&]() {
		field_.overlay_par(lphs_factors, multiply<Complex, Complex>);
		};



	for (unsigned step = 0; step < profile_.frame_interval; ++step) {

		//TSSP substep1
		tssp1();

		//TSSP substep2
		fft2d_ortho_par(field_);
		tssp2();
		ifft2d_ortho_par(field_);

		//TSSP substep3
		tssp1();
	}


	//-------------------------------------render--------------------------------------

	CompressedHeatMap heatmap(0.6, 1.0, Color::Sandstone());
	const auto colorizer = [&](const Complex& c) {
		return heatmap(c.amplitude());
		};

	mat<RGB> field_img = mat<RGB>::creat_par(field_, colorizer);


	++Sim::step_cur_;
	return field_img;
}

CGLESim::~CGLESim()
{
	writeMat(WF_FILE_PATH, field_, profile_);
}