#include "plot_test.h"

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
#include <numer_cube.h>
#include <numer_file.h>
#include <functional>
#include <vector>
#include <algorithm>

#include <stbi.h>


using namespace numer;





void PlotTest::run()
{

	std::cout << "\n>..PlotTest: begin" << std::endl;
	BENCHMARK_BEGIN(plot_test);


	HydrogenState nlm(4, 3, 0, 0.02);
	Cartes3Adp nlm_in_xyz(nlm);
	

	BENCHMARK_BEGIN(render);

	mat<RGB> wf3d_img = DensityPlot3D(720, 1280)
		.setCamAzimuthAngle(0.4 * Pi)
		.setCamElevationAngle(0.1 * Pi)
		.setBrightnessGain(3.0)
		.setFineness(500)
		.renderImage(nlm_in_xyz, std::function<Vec3<double>(Complex)>(gridclr::HeatMap(0.0, 2.0, LinearHeatMap(0.0, 2.0, Color::Zone_odd()))));

	BENCHMARK_END(render);



	stbi::ImageWriter<stbi::format::PNG> writer;
	
	writer.writeInto("./image/plot/wavefunc", wf3d_img[0], wf3d_img.ncols(), wf3d_img.nrows());




	std::cout << "\n>..PlotTest: end" << std::endl;
	BENCHMARK_END(plot_test);
}