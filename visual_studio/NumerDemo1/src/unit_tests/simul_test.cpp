#include "simul_test.h"

#include <string>
#include <iostream>
#include "benchmark.h"
#include "simul_utils.h"
#include "processbar.h"
#include <numer_mat.h>
#include <numer_visualize.h>
#include <stbi.h>

#include <memory>
#include "sim/sim.h"
#include "sim/ho1d_sim.h"

using namespace numer;



constexpr unsigned FRAMES = 180;
constexpr unsigned START_FRAME = 00;




void SimulTest::run()
{
	std::cout << "\n>..SimulTest: begin" << std::endl;
	BENCHMARK_BEGIN(simuul_test);

	SimImagePath path_keeper("ho1d", START_FRAME);

	BENCHMARK_BEGIN(pre_computation);
	std::unique_ptr<Sim> up_sim = std::make_unique<HO1dSim>(FRAMES);
	BENCHMARK_END(pre_computation);

	stbi::ImageWriter<stbi::format::PNG> img_w;


	std::atomic<unsigned int> counter(0);
	std::atomic<float> cycle_timer(0.0);
	Tracker tk_process("Simulation process:", counter, cycle_timer, FRAMES);

	for (unsigned i = 0; i < FRAMES; ++i) {
		auto start = std::chrono::high_resolution_clock::now();

		mat<RGB> frame = up_sim->renderFrame();
		img_w.writeInto(path_keeper.getPath(i), frame[0], frame.ncols(), frame.nrows());

		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		float dura_sec = (float)duration / 1000.0f;
		counter.fetch_add(1, std::memory_order_relaxed);
		cycle_timer.store(dura_sec, std::memory_order_relaxed);
	}

	tk_process.stop();

	std::cout << "\n>..SimulTest: end" << std::endl;
	BENCHMARK_END(simuul_test);
}