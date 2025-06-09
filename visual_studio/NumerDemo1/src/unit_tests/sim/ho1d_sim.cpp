#include "ho1d_sim.h"

#include <numer_complex.h>
#include <numer_visualize.h>
#include <numer_mat.h>
#include <numer_fourier.h>
#include <numer_eigenfunc.h>
#include <numer_grid.h>
#include <numer_plot.h>
#include <numer_indexer.h>
#include <functional>
#include <algorithm>
#include <execution>





using namespace numer;



constexpr size_t N = 3000;
constexpr double L = 20.0;
constexpr RangeSpec SPEC{ -L, L, N };
constexpr size_t M = 151;
constexpr double omega_t_total = 2.0 * Pi;


constexpr size_t HEIGHT = 2160;
constexpr size_t WIDTH = 3840;

constexpr size_t hght_u = 200;
constexpr size_t wdth = 3000;
constexpr RGB bg_clr{ 0, 0, 0 };
constexpr RGB line_clr{ 225, 225, 225 };




HO1dSim::HO1dSim(unsigned Steps)
	: Sim(Steps), nbasis_bra_(M), nket_(M), q_distri_gen_(100, 64)
{
	for (size_t i = 0; i < M; ++i) {
		std::vector<Complex> phi_i(N);
		UnaryFuncSampler<HermiFunc> sampler(HermiFunc(i, 1.0), SPEC);
		for (size_t j = 0; j < N; ++j) {
			phi_i[j] = sampler(j);
		}
		fft_ortho(phi_i, N);
		std::for_each(phi_i.begin(), phi_i.end(), [](Complex& C) {C = C.conj(); });
		nbasis_bra_[i] = std::move(phi_i);
	}



	CoherentState1D phi_a(2.5 * Complex::expi(0.0 * Pi), 0.5);

	auto oddCohrent = [](double x)->Complex {
		static CoherentState1D phia(20.0 * Complex::expi(0.0 * Pi), 4.0);
		static CoherentState1D phib(1.25 * Complex::expi(0.0 * Pi), 0.25);
		return (phia(x) + phib(x)) / sqrt(2.0);
		};


	UnaryFuncSampler<std::function<Complex(double)>> samp1(oddCohrent, SPEC);
	std::vector<Complex> ket(N);
	for (unsigned i = 0; i < N; ++i) {
		ket[i] = samp1(i);
	}
	fft_ortho(ket, N);
	for (size_t i = 0; i < M; ++i) {
		nket_[i] = std::inner_product(ket.cbegin(), ket.cend(), nbasis_bra_[i].cbegin(), Complex::zero()) * samp1.diffElem();
	}

}





mat<RGB> HO1dSim::renderFrame() {

	//----------------------------------time parameter---------------------------------

	double t_para = static_cast<double>(Sim::step_cur_) / static_cast<double>(Sim::steps_);
	double omega_t = omega_t_total * t_para;
	//double beta_t = 1.0 + 4.0 * t_para;
	
	//-------------------------------------simulation----------------------------------

/*
	CoherentState1D phi_a(0.0 * Complex::expi(0.0 * Pi), beta_t);

	UnaryFuncSampler<std::function<Complex(double)>> samp1(phi_a, SPEC);
	std::vector<Complex> rket1(N);
	for (unsigned i = 0; i < N; ++i) {
		rket1[i] = samp1(i);
	}
	std::vector<Complex> pket1 = rket1, nket1(M);
	fft_ortho(pket1, N);
	for (size_t i = 0; i < M; ++i) {
		nket1[i] = std::inner_product(pket1.cbegin(), pket1.cend(), nbasis_bra_[i].cbegin(), Complex::zero()) * samp1.diffElem();
	}
*/






	std::vector<Complex> nket1 = nket_;

	Complex evo = Complex::expi(omega_t);
	Complex phashft_n = Complex::expi(omega_t / 2.0);
	for (size_t i = 0; i < M; ++i) {
		nket1[i] *= phashft_n;
		phashft_n *= evo;
	}

	std::vector<Complex> rket1(N, Complex::zero());
	std::vector<size_t> args(N);
	for (size_t i = 0; i < N; ++i) {
		args[i] = i;
	}
	std::for_each(std::execution::par, args.begin(), args.end(), [&](size_t id) {
		for (size_t j = 0; j < M; ++j) {
			rket1[id] += nket1[j] * nbasis_bra_[j][id].conj();
		}
	});

	ifft_ortho(rket1, N);


	//-------------------------------------render-----------------------------------


	mat<RGB> plot_img(HEIGHT, WIDTH);


	RangeSampler y_idxer1(RangeSpec{ -1.0, 1.0, hght_u * 4 });
	RangeSampler y_idxer2(RangeSpec{ -0.1, 0.4, hght_u * 4});
	RangeSampler x_idxer(RangeSpec{ 0, N - 1, wdth });
	RangeSampler n_idxer(RangeSpec{ 0, M - 1, wdth });

	constexpr size_t y_off0 = 80ULL + 0;
	constexpr size_t y_off1 = 80ULL + hght_u * 4;
	constexpr size_t y_off2 = 80ULL + hght_u * 5;
	constexpr size_t y_off3 = 80ULL + hght_u * 9;


	ComplxRainbowClr rvs(0.5);
	ComplxRainbowClr nvs(0.15);

	mat<RGB> rspec = HeatMapPlot(hght_u, wdth).renderImage(VecIndexer(rket1), rket1.size(), std::function<RGB(Complex)>(rvs));
	mat<RGB> nspec = HeatMapPlot(50, wdth - 1100).renderImage(VecIndexer(nket1), 101, std::function<RGB(Complex)>(nvs));




	LinearHeatMap heatmap1(0.0, 0.25, Color::Viridis());
	LinearHeatMap heatmap2(-0.1, 1.0, Color::Vaporwave());

	auto ComplxAmpClr1 = [&heatmap1](Complex C) -> RGB {
		return heatmap1(C.sqrdAmp());
		};

	auto ComplxAmpClr2 = [&heatmap2](Complex C) -> RGB {
		return heatmap2(C.sqrdAmp());
		};

	mat<RGB> rket_plot = Histogram(hght_u * 4, wdth, -1.2, 1.2)
		.setBackColor(RGB{ 0, 0, 0 })
		.setLineTrans(127)
		.setLineColor(RGB{ 0xf1, 0xbe, 0x5e })
		.drawData(VecIndexer(rket1), rket1.size(), [](Complex& C) {return C.im(); }, std::function<RGB(Complex)>(ComplxAmpClr2))
		.drawData(VecIndexer(rket1), rket1.size(), [](Complex& C) {return C.re(); }, std::function<RGB(Complex)>(ComplxAmpClr1))
		.setLineColor(RGB{ 225, 225, 225 })
		.drawHorizLine(0.0)
		.getImage();

	mat<RGB> nket_plot = Histogram(hght_u * 4, wdth - 1100, -0.1, 0.4)
		.setBackColor(RGB{ 0, 0, 0 })
		.drawData(VecIndexer(nket1), 101, [](Complex& C) {return C.sqrdAmp(); })
		.drawHorizLine(0.0)
		.getImage();




	mat<double> q_distri = q_distri_gen_.calculate(VecIndexer(nket1), nket1.size());
	MatIndexer q_distri_v(q_distri);

	mat<RGB> hq_img = HeatMapPlot(800, 800).renderImage(ImageOrientation(q_distri_v, q_distri.ncols(), q_distri.nrows()), q_distri.nrows(), q_distri.ncols(), std::function<RGB(double)>(LinearHeatMap(-0.01, 0.1, Color::Vaporwave())));


	plot_img.overlay(rket_plot, y_off0, 420, [](RGB cl1, RGB cl2) {return cl2; });
	plot_img.overlay(rspec, y_off1, 420, [](RGB cl1, RGB cl2) {return cl2; });
	plot_img.overlay(nket_plot, y_off2 + 50, 420, [](RGB cl1, RGB cl2) {return cl2; });
	plot_img.overlay(hq_img, y_off2 + 100, 420 + 2200, [](RGB cl1, RGB cl2) {return cl2; });
	plot_img.overlay(nspec, y_off3 + 50, 420, [](RGB cl1, RGB cl2) {return cl2; });

	++Sim::step_cur_;
	return plot_img;
}
