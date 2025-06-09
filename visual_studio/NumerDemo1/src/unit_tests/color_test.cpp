#include "color_test.h"

#include "benchmark.h"
#include <iostream>
#include <string>
#include <numer_visualize.h>
#include <numer_mat.h>
#include <stbi.h>

using namespace numer;






static inline
numer::RGB linearMix(double ratio_1, const numer::RGB& first, const numer::RGB& second)
{
	return numer::RGB{
		(uint8_t)(ratio_1 * first.R + (1.0 - ratio_1) * second.R),
		(uint8_t)(ratio_1 * first.G + (1.0 - ratio_1) * second.G),
		(uint8_t)(ratio_1 * first.B + (1.0 - ratio_1) * second.B)
	};
}

static numer::RGB clrNew(double X) {

	constexpr double t0 = 0.0, t1 = 0.4, t2 = 0.55, t3 = 0.7, t4 = 0.8, t5 = 0.9, t6 = 1.0;
	constexpr numer::RGB
		Clr0{ 0x00, 0x00, 0x00 },
		Clr1{ 0x71, 0x21, 0x0c },
		Clr2{ 0xb3, 0x3b, 0x1b },
		Clr3{ 0xe6, 0x98, 0x42 },
		Clr4{ 0xe6, 0xd4, 0x65 },
		Clr5{ 0xea, 0xe6, 0xbc },
		Clr6{ 0xf1, 0xfd, 0xff };


	if (X < t0) return Clr0;
	else if (X < t1) return linearMix((t1 - X) / (t1 - t0), Clr0, Clr1);
	else if (X < t2) return linearMix((t2 - X) / (t2 - t1), Clr1, Clr2);
	else if (X < t3) return linearMix((t3 - X) / (t3 - t2), Clr2, Clr3);
	else if (X < t4) return linearMix((t4 - X) / (t4 - t3), Clr3, Clr4);
	else if (X < t5) return linearMix((t5 - X) / (t5 - t4), Clr4, Clr5);
	else if (X < t6) return linearMix((t6 - X) / (t6 - t5), Clr5, Clr6);
	else return Clr6;
}


void ColorTest::run()
{
	stbi::ImageWriter<stbi::format::PNG> writer;
	stbi::ImageLoader loader;


	std::cout << "\n>..ColorTest: begin" << std::endl;
	BENCHMARK_BEGIN(color_test);

	if (loader.load("./image/color/original/gradient.png", stbi::LOAD_GRY))
	{
		mat<uint8_t> img(loader.height(), loader.width(), 0);
		loader.putInto(img.begin());

		writer.writeInto("./image/color/grayscale", img.cbegin(), img.ncols(), img.nrows());

		auto img_num = mat<double>::creat_par(img, GrayScale(0.0, 1.0));

		auto img_re = mat<RGB>::creat_par(img_num, Color::Rainbow());

		writer.writeInto("./image/color/newcolor", img_re[0], img_re.ncols(), img_re.nrows());
	}
	
	std::cout << "\n>..ColorTest: end" << std::endl;
	BENCHMARK_END(color_test);
}
