#include "fourier3_test.h"

#include "benchmark.h"
#include <iostream>
#include <string>
#include <numer_complex.h>
#include <numer_matrix.h>
#include <numer_visualize.h>
#include <numer_mat.h>
#include <numer_fourier.h>
#include <stbi.h>//requires stbi library: https://github.com/modulelife/stbi




using namespace numer;



void Fourier3Test::run()
{

	std::cout << "\n>..Fourier3Test: begin" << std::endl;
	BENCHMARK_BEGIN(fourier3_test);



	stbi::ImageWriter<stbi::format::PNG> writer;
	stbi::ImageLoader loader;


	loader.load("./image/fourier/testimg/chongqing.jpg", stbi::LOAD_RGB);
	mat<RGB> img(loader.height(), loader.width());
	loader.putInto(img.begin());

	writer.writeInto("./image/fourier/original", img.cbegin(), img.ncols(), img.nrows());

	
	using cvec3 = Vec3<Complex>;
	GrayScale double_byte(0.0, 1.0);

	const auto RGB_to_cvec3 = [&double_byte](const RGB& pix) {
		cvec3 vpix{};
		vpix[0] = double_byte(pix.R);
		vpix[1] = double_byte(pix.G);
		vpix[2] = double_byte(pix.B);
		return vpix;
		};

	const auto cvec3_to_RGB = [&double_byte](const cvec3& vpix) {
		RGB pix{};
		pix.R = double_byte(vpix[0].amplitude());
		pix.G = double_byte(vpix[1].amplitude());
		pix.B = double_byte(vpix[2].amplitude());
		return pix;
		};

	auto vimg = mat<cvec3>::creat(img, RGB_to_cvec3);
	fft2d_ortho_par(vimg);

	auto amp_img = mat<RGB>::creat(vimg, cvec3_to_RGB);
	centralize(amp_img);

	ifft2d_ortho_par(vimg);

	auto img_re = mat<RGB>::creat(vimg, cvec3_to_RGB);

	writer.writeInto("./image/fourier/amp_spectr3", amp_img.cbegin(), amp_img.ncols(), amp_img.nrows());
	writer.writeInto("./image/fourier/reconstructed3", img_re.cbegin(), img_re.ncols(), img_re.nrows());




	std::cout << "\n>..Fourier3Test: end" << std::endl;
	BENCHMARK_END(fourier3_test);
}