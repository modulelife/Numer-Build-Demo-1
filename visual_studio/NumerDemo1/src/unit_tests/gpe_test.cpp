#include "gpe_test.h"

#include "benchmark.h"
#include <iostream>
#include <string>

#include <numer_visualize.h>
#include <numer_mat.h>
#include <numer_fourier.h>
#include <numer_eigenfunc.h>
#include <numer_grid.h>
#include <numer_plot.h>
#include <numer_qmkit.h>
#include <numer_indexer.h>
#include <numer_matrix.h>
#include <numer_cube.h>
#include <numer_file.h>
#include <functional>
#include <vector>
#include <algorithm>

#include <stbi.h>


using namespace numer;



constexpr unsigned N = 256;
constexpr double L = 20.0;
constexpr RangeSpec spec{ -L, L, N };

constexpr Complex alpha = 1.5 * Complex::identity();
constexpr double gy = 1.0, gz = 1.0;
static CoherentState3D trial_wf(0.0, 1.0, 0.0, gy, 0.0, gy);

static const CoherentState3D drop1(0.0, 1.0, 0.0, gy, -1.5, gz);
static const CoherentState3D drop2(0.0, 1.0, 0.5 * Complex::i(), gy, 1.5, gz);
static const auto init_wf = [](double x, double y, double z) {
	return drop1(x, y, z) + drop2(x, y, z);
	};



constexpr double E = 100000.0, d = 2.0, s = 1.8;


void GPETest::run()
{

	std::cout << "\n>..GPETest: begin" << std::endl;
	BENCHMARK_BEGIN(gpe_test);



	const auto vint_func = [=](double x, double y, double z) {
		constexpr double V0 = 10.0, R = 5.0;
		double r = sqrt(x * x + y * y + z * z);
		if (r == 0.0) return 0.0;
		//double vint = 4.0 * E * (pow(d / (r + s), 12) - pow(d / (r + s), 6));
		//return Complex{ vint, 0.0 };
		double r_xy = sqrt(x * x + y * y);
		double phi = atan2(y, x); // 方位角
		double twist = cos(5 * phi); // 五重对称性
		return V0 * exp(-(r_xy * r_xy + z * z) / (R * R)) * twist;
		};

	const TernaryFuncSampler vint_samp(vint_func, spec, spec, spec);

	auto vintsp = cube<Complex>::creat_par(N, N, N, vint_samp);
	//fft3d_ortho_par(vintsp);
	//centralize(vintsp);

	FieldRelocator field(vintsp, N, N, N);
	
	mat<RGB> vintsp3d_img = DensityPlot3D(500, 500)
		.setCamAzimuthAngle(0.1 * Pi)
		.setCamElevationAngle(0.05 * Pi)
		.setBrightnessGain(1.0)
		.setFineness(500)
		.renderImage(field, std::function<Vec3<double>(Complex)>(gridclr::HeatMap(0.0, 1.0, Color::Glacier())));



	std::cout << vint_func(0.0, 0.0, 0.0) << std::endl;


	stbi::ImageWriter<stbi::format::PNG> writer;

	writer.writeInto("./image/gpe/vintsp", vintsp3d_img[0], vintsp3d_img.ncols(), vintsp3d_img.nrows());






	std::cout << "\n>..GPETest: end" << std::endl;
	BENCHMARK_END(gpe_test);
}
