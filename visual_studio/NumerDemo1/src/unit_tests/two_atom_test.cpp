#include "two_atom_test.h"

#include "benchmark.h"
#include <iostream>
#include <string>
#include <numer_complex.h>
#include <numer_visualize.h>
#include <numer_mat.h>
#include <numer_fourier.h>
#include <numer_eigenfunc.h>
#include <numer_grid.h>
#include <numer_plot.h>
#include <numer_qmkit.h>
#include <numer_indexer.h>
#include <numer_matrix.h>
#include <functional>
#include <vector>
#include <algorithm>

#include <stbi.h>


using namespace numer;





const auto J0 = [](double r) -> double {
    return std::cyl_bessel_j(0, r);
    };

const auto J1 = [](double r) -> double {
    return std::cyl_bessel_j(1, r);
    };

const auto Jc1 = [](double r) -> double {
    return std::cyl_bessel_j(1, r) / r;
    };

template<class Func>
static double integral(const RangeSpec& Range, Func&& Function) {
	std::vector<double> sampls(Range.length);
	UnaryFuncSampler samplr(Function, Range);
	for (unsigned i = 0; i < Range.length; ++i) {
		sampls[i] = samplr(i);
	}
	return samplr.diffElem() * std::accumulate(sampls.cbegin(), sampls.cend(), 0.0);
}

template<typename Ty, unsigned Dim>
constexpr auto approx_expm(const SquareMatrix<Ty, Dim>& M_) {
	constexpr unsigned max_n = 10;
	SquareMatrix<Ty, Dim> I = makeIdentityMatrix<Ty, Dim>();
	SquareMatrix<Ty, Dim> result = I;
	for (unsigned i = max_n; i > 0; --i) {
		result = I + M_ * result / static_cast<double>(i);
	}
	return result;
}




constexpr Vec3<double> ex{ 1.0, 0.0, 0.0 }, ey{ 0.0, 1.0, 0.0 }, ez{ 0.0, 0.0, 1.0 };
constexpr auto I = makeIdentityMatrix<double, 3>();
constexpr auto xx = ex * ex.t(), yy = ey * ey.t(), zz = ez * ez.t(), zx = ez * ex.t(), xz = ex * ez.t();

constexpr double z1 = 0.2, z2 = 0.2, x12 = 3.0;
constexpr Vec3<double> r1{ 0.0, 0.0, z1 }, r2{ x12, 0.0, z2 }, r12 = r1 - r2;
constexpr Vec3<double> d1{ 1.0, 1.0, 1.0 }, d2{ 1.0, 1.0, 1.0 };
constexpr RangeSpec theta{ 0.001, Pi / 2.0, 2000 };

constexpr double TIME = 5.0;
constexpr unsigned STEPS = 2000;
constexpr double dt = TIME / STEPS;



void TwoAtomTest::run()
{

	std::cout << "\n>..TwoAtomTest: begin" << std::endl;
	BENCHMARK_BEGIN(two_atom_test);



	BENCHMARK_BEGIN(matrix_element);
	
	const Vec3<double> ed1 = d1 / sqrt(d1.t() * d1);
	const Vec3<double> ed2 = d2 / sqrt(d2.t() * d2);
	const double r_12 = sqrt(r12.t() * r12);
	const Vec3<double> er12 = r12 / r_12;
	const auto rr = er12 * er12.t();



	const double gs_1 =
		0.75 * ed1.t() * (-I * (sin(2 * z1) / (2 * z1) + cos(2 * z1) / pow(2 * z1, 2) - sin(2 * z1) / pow(2 * z1, 3)) + zz * (sin(2 * z1) / (2 * z1) - cos(2 * z1) / pow(2 * z1, 2) + sin(2 * z1) / pow(2 * z1, 3))) * ed1;

	const double gs_2 =
		0.75 * ed2.t() * (-I * (sin(2 * z2) / (2 * z2) + cos(2 * z2) / pow(2 * z2, 2) - sin(2 * z2) / pow(2 * z2, 3)) + zz * (sin(2 * z2) / (2 * z2) - cos(2 * z2) / pow(2 * z2, 2) + sin(2 * z2) / pow(2 * z2, 3))) * ed2;

	const double ds_1 =
		0.75 * ed1.t() * (-I * (cos(2 * z1) / (2 * z1) + sin(2 * z1) / pow(2 * z1, 2) + cos(2 * z1) / pow(2 * z1, 3)) + zz * (cos(2 * z1) / (2 * z1) + 3.0 * sin(2 * z1) / pow(2 * z1, 2) + 3.0 * cos(2 * z1) / pow(2 * z1, 3))) * ed1;

	const double ds_2 =
		0.75 * ed2.t() * (-I * (cos(2 * z2) / (2 * z2) + sin(2 * z2) / pow(2 * z2, 2) + cos(2 * z2) / pow(2 * z2, 3)) + zz * (cos(2 * z2) / (2 * z2) + 3.0 * sin(2 * z2) / pow(2 * z2, 2) + 3.0 * cos(2 * z2) / pow(2 * z2, 3))) * ed2;

	const Complex A11 = -(0.5 + gs_1) + Complex::i() * ds_1;
	const Complex A22 = -(0.5 + gs_2) + Complex::i() * ds_2;

	

	const double gcd =
		0.75 * ed1.t() * (I * (sin(r_12) / r_12 + cos(r_12) / pow(r_12, 2) - sin(r_12) / pow(r_12, 3)) - rr * (sin(r_12) / r_12 + 3.0 * cos(r_12) / pow(r_12, 2) - 3.0 * sin(r_12) / pow(r_12, 3))) * ed2;

	const double dcd =
		0.75 * ed1.t() * (I * (cos(r_12) / r_12 - sin(r_12) / pow(r_12, 2) - cos(r_12) / pow(r_12, 3)) - rr * (cos(r_12) / r_12 - 3.0 * sin(r_12) / pow(r_12, 2) - 3.0 * cos(r_12) / pow(r_12, 3))) * ed2;


	const double gcr = 0.75 * ed1.t() * (
		-(xx - yy) * integral(theta, [](double t) {return pow(sin(t), 3) * Jc1(x12 * sin(t)) * cos((z1 + z2) * cos(t)); })
		- (yy - zz) * integral(theta, [](double t) {return sin(t) * J0(x12 * sin(t)) * cos((z1 + z2) * cos(t)); })
		- (xx + zz) * integral(theta, [](double t) {return sin(t) * pow(cos(t), 2) * J0(x12 * sin(t)) * cos((z1 + z2) * cos(t)); })
		+ (zx - xz) * integral(theta, [](double t) {return pow(sin(t), 2) * cos(t) * J1(x12 * sin(t)) * sin((z1 + z2) * cos(t)); })) * ed2;

	const auto FJ =
		-(xx - yy) * (Jc1(x12) / (z1 + z2) + 2.0 * Jc1(x12) / pow(z1 + z2, 3))
		- (yy - zz) * (J0(x12) / (z1 + z2))
		+ (xx + zz) * (2.0 * J0(x12) / pow(z1 + z2, 3))
		+ (zx - xz) * (J1(x12) / pow(z1 + z2, 2) + 6.0 * J1(x12) / pow(z1 + z2, 4));


	const double dcr = 0.75 * ed1.t() * (
		(xx - yy) * integral(theta, [](double t) {return pow(sin(t), 3) * Jc1(x12 * sin(t)) * sin((z1 + z2) * cos(t)); })
		+ (yy - zz) * integral(theta, [](double t) {return sin(t) * J0(x12 * sin(t)) * sin((z1 + z2) * cos(t)); })
		+ (xx + zz) * integral(theta, [](double t) {return sin(t) * pow(cos(t), 2) * J0(x12 * sin(t)) * sin((z1 + z2) * cos(t)); })
		+ (zx - xz) * integral(theta, [](double t) {return pow(sin(t), 2) * cos(t) * J1(x12 * sin(t)) * cos((z1 + z2) * cos(t)); })
		+ FJ) * ed2;

	const Complex A12 = -(gcd + gcr) + Complex::i() * (dcd + dcr);
	const Complex A21 = A12;


	std::cout << gs_1 << std::endl;
	std::cout << gs_2 << std::endl;
	std::cout << ds_1 << std::endl;
	std::cout << ds_2 << std::endl;
	std::cout << gcd << std::endl;
	std::cout << dcd << std::endl;
	std::cout << gcr << std::endl;
	std::cout << dcr << std::endl;



	BENCHMARK_END(matrix_element);




	BENCHMARK_BEGIN(evolution);

	Carre<Complex, 2> A{
		std::array{A11, A12},
		std::array{A21, A22}
	};

	Carre<Complex, 2> U = expm_approx<10>(A * dt);
	Vec<Complex, 2> C{ 1.0, 0.0 };

	std::vector<Vec<Complex, 2>> evo(STEPS);

	for (unsigned i = 0; i < STEPS; ++i) {
		evo[i] = C;
		C = U * C;
	}


	BENCHMARK_END(evolution);




	BENCHMARK_BEGIN(ploting);

	std::vector<double> prob1(STEPS);
	std::vector<double> prob2(STEPS);
	std::vector<double> eri(STEPS);
	std::vector<double> ercpi1(STEPS);
	std::vector<double> ercpi2(STEPS);
	std::vector<double> probs(STEPS);
	std::vector<double> proba(STEPS);
	std::vector<double> probdiff(STEPS);

	for (unsigned i = 0; i < STEPS; ++i) {
		prob1[i] = evo[i][0].sqrdAmp();
		prob2[i] = evo[i][1].sqrdAmp();
		eri[i] = 2 * Re(Complex::i() * (A12 * evo[i][0].conj() * evo[i][1] + A21 * evo[i][0] * evo[i][1].conj()));
		ercpi1[i] = -2 * A11.im() * evo[i][0].sqrdAmp();
		ercpi2[i] = -2 * A22.im() * evo[i][1].sqrdAmp();
		probs[i] = (evo[i][0] + evo[i][1]).sqrdAmp() / 2.0;
		proba[i] = (evo[i][0] - evo[i][1]).sqrdAmp() / 2.0;
		probdiff[i] = probs[i] - proba[i];
	}


	mat<RGB> prob_img = Histogram(1000, 1500, 0.0, 1.0)
		.setBackColor(RGB{ 225, 225, 225 })
		.setLineColor(RGB{ 0, 0, 0 })
		.drawHorizLine(0.0)
		.setLineTrans(127)
		.setLineColor(RGB{ 125, 177, 251 })
		.drawData(VecIndexer(prob1), STEPS)
		.setLineColor(RGB{ 253, 120, 110 })
		.drawData(VecIndexer(prob2), STEPS)
		.getImage();

	mat<RGB> e_img = Histogram(1000, 1500, 0.0, 2.5)
		.setBackColor(RGB{ 225, 225, 225 })
		.setLineColor(RGB{ 0, 0, 0 })
		.drawHorizLine(0.0)
		.setLineTrans(127)
		.setLineColor(RGB{ 253, 70, 70 })
		.drawData(VecIndexer(ercpi1), STEPS)
		.setLineColor(RGB{ 50, 160, 50 })
		.drawData(VecIndexer(ercpi2), STEPS)
		.setLineColor(RGB{ 10, 10, 170 })
		.drawData(VecIndexer(eri), STEPS)
		.getImage();

	mat<RGB> psa_img = Histogram(1000, 1500, -0.1, 0.5)
		.setBackColor(RGB{ 225, 225, 225 })
		.setLineColor(RGB{ 0, 0, 0 })
		.drawHorizLine(0.0)
		.setLineTrans(127)
		.setLineColor(RGB{ 253, 70, 70 })
		.drawDataLine(VecIndexer(proba), STEPS)
		.setLineColor(RGB{ 10, 10, 170 })
		.drawDataLine(VecIndexer(probs), STEPS)
		.setLineColor(RGB{ 50, 160, 50 })
		.drawDataLine(VecIndexer(probdiff), STEPS)
		.getImage();


	BENCHMARK_END(ploting);


	stbi::ImageWriter<stbi::format::PNG> writer;


	writer.writeInto("./image/two_atom/probability", prob_img[0], prob_img.ncols(), prob_img.nrows());
	writer.writeInto("./image/two_atom/energy", e_img[0], e_img.ncols(), e_img.nrows());
	writer.writeInto("./image/two_atom/psa", psa_img[0], psa_img.ncols(), psa_img.nrows());


	std::cout << "\n>..TwoAtomTest: end" << std::endl;
	BENCHMARK_END(two_atom_test);
}