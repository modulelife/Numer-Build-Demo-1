#include <iostream>
#include <memory>
#include <vector>
#include "unit_tests/test.h"
#include "unit_tests/color_test.h"
#include "unit_tests/plot_test.h"
#include "unit_tests/simul_test.h"
#include "unit_tests/mandelbrot_test.h"
#include "unit_tests/spin_test.h"
#include "unit_tests/fourier3_test.h"






int main()
{
	std::unique_ptr<Test> up_test;
	unsigned op;
	std::cout << "\n\tchoose a test:\n\t\t1.ColorTest\n\t\t2.PlotTest\n\t\t3.SimulTest\n\t\t4.MandelbrotTest\n\t\t5.SpinTest\n\t\t6.Fourier3Test\n\n\tenter: ";
	std::cin >> op;

	switch (op) {
	case 1:
		up_test = std::make_unique<ColorTest>();
		break;
	case 2:
		up_test = std::make_unique<PlotTest>();
		break;
	case 3:
		up_test = std::make_unique<SimulTest>();
		break;
	case 4:
		up_test = std::make_unique<MandelbrotTest>();
		break;
	case 5:
		up_test = std::make_unique<SpinTest>();
		break;
	case 6:
		up_test = std::make_unique<Fourier3Test>();
		break;
	default: 
		up_test = std::make_unique<Test>();
		break;
	}

	up_test->run();

	return 0;
}