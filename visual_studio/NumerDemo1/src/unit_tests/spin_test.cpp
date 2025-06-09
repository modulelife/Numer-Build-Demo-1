#include "spin_test.h"

#include "benchmark.h"
#include <iostream>
#include <string>
#include <numer_complex.h>
#include <numer_mat.h>
#include <numer_plot.h>
#include <numer_indexer.h>
#include <numer_matrix.h>
#include <numer_fourier.h>
#include <vector>

#include <stbi.h>


using namespace numer;


constexpr Carre<Complex, 2> pauli_x{
	std::array{Complex::zero(), Complex::identity()},
	std::array{Complex::identity(), Complex::zero()}
};

constexpr Carre<Complex, 2> pauli_y{
	std::array{Complex::zero(), -Complex::i()},
	std::array{Complex::i(), Complex::zero()}
};

constexpr Carre<Complex, 2> pauli_z{
	std::array{Complex::identity(), Complex::zero()},
	std::array{Complex::zero(), -Complex::identity()}
};

constexpr Vec3<Carre<Complex, 2>> pauli_vec{ pauli_x, pauli_y, pauli_z };

constexpr double B_0 = 0.5;
constexpr double B_1 = 1.0;

constexpr double omega = 0.1;



constexpr double TIME = 5000.0;
constexpr unsigned STEPS = 10000;
constexpr double dt = TIME / STEPS;

const auto t_switch = [](double time) {
    if (time > 100.0 && time < 4000.0) return 1.0;
    else return 0.0;
    };



void SpinTest::run()
{

	std::cout << "\n>..SpinTest: begin" << std::endl;
	BENCHMARK_BEGIN(spin_test);

    
    Vec<Complex, 2> spinor{
        Complex::identity(),
        Complex::zero()
    };

    std::vector<double> time_points(STEPS);
    std::vector<double> s_x(STEPS), s_y(STEPS), s_z(STEPS);
    std::vector<Complex> s_tran(STEPS);
    std::vector<double> prob_up(STEPS), prob_down(STEPS);

    for (unsigned step = 0; step < STEPS; ++step) {

        const double t = step * dt;
        Vec3<double> B_vec{ 
            B_1 * cos(omega * t),
            B_1 * sin(omega * t),
            B_0
        };
        //hamiltonian H(t) = -1/2 B(t)·σ
        auto H = (-0.5) * B_vec.t() * pauli_vec;
        //evolution operator U = exp(-iHdt)
        auto U = expm_approx<10>((-Complex::i()) * H * dt);
        //evolve the ket ψ(t+dt) = U(t+dt,t)ψ(t)
        spinor = U * spinor;
        // ⟨σ⟩ = ψ†σψ
        double sx = Re(spinor.hc() * pauli_x * spinor);
        double sy = Re(spinor.hc() * pauli_y * spinor);
        s_x[step] = sx;
        s_y[step] = sy;
        s_tran[step] = Complex{ sx, sy };
        s_z[step] = Re(spinor.hc() * pauli_z * spinor);
        // probability |⟨↑|ψ⟩|² & |⟨↓|ψ⟩|²
        prob_up[step] = spinor[0].sqrdAmp();
        prob_down[step] = spinor[1].sqrdAmp();
        time_points[step] = t;
    }

    mat<RGB> spin_img(1000, STEPS);

    mat<RGB> prob_plot = Histogram(200, STEPS, -0.1, 1.1)
        .setBackColor(RGB{ 225, 225, 225 })
        .setLineColor(RGB{ 0, 0, 0 })
        .drawHorizLine(0.0)
        .setLineTrans(127)
        .setLineColor(RGB{ 253, 120, 110 })
        .drawData(VecIndexer(prob_up), STEPS)
        .setLineColor(RGB{ 125, 177, 251 })
        .drawData(VecIndexer(prob_down), STEPS)
        .getImage();
    mat<RGB> sxy_bar = HeatMapPlot(50, STEPS)
        .renderImage(VecIndexer(s_tran), STEPS, std::function<RGB(Complex)>(ComplxPhaseClr(Color::Glacier())));
    mat<RGB> sxy_plot = Histogram(300, STEPS, -1.2, 1.2)
        .setBackColor(RGB{ 225, 225, 225 })
        .setLineColor(RGB{ 0, 0, 0 })
        .drawHorizLine(0.0)
        .setLineTrans(127)
        .setLineColor(RGB{ 253, 70, 70 })
        .drawDataLine(VecIndexer(s_x), STEPS)
        .setLineColor(RGB{ 10, 10, 170 })
        .drawDataLine(VecIndexer(s_y), STEPS)
        .getImage();
    mat<RGB> sz_plot = Histogram(300, STEPS, -1.2, 1.2)
        .setBackColor(RGB{ 225, 225, 225 })
        .setLineColor(RGB{ 0, 0, 0 })
        .drawHorizLine(0.0)
        .setLineTrans(127)
        .setLineColor(RGB{ 50, 160, 50 })
        .drawDataLine(VecIndexer(s_z), STEPS)
        .getImage();
    mat<RGB> sz_bar = HeatMapPlot(50, STEPS)
        .renderImage(VecIndexer(s_z), STEPS, std::function<RGB(double)>(LinearHeatMap(-1.0, 1.0, Color::Zone_half())));

    fft_ortho(s_tran, STEPS);
    LinearHeatMap psd_heatmap(0.0, 10.0, Color::Thermo());
    const auto psd_clr = [&psd_heatmap](Complex c) {return psd_heatmap(c.sqrdAmp()); };
    mat<RGB> psd_bar = HeatMapPlot(100, STEPS)
        .renderImage(VecReIndexer(s_tran, id_order::cyclic(STEPS, STEPS / 2)), STEPS, std::function(psd_clr));


    spin_img.overlay(prob_plot, 0, 0, chooseSecond<RGB>)
        .overlay(psd_bar, 200, 0, chooseSecond<RGB>)
        .overlay(sxy_bar, 300, 0, chooseSecond<RGB>)
        .overlay(sxy_plot, 350, 0, chooseSecond<RGB>)
        .overlay(sz_bar, 650, 0, chooseSecond<RGB>)
        .overlay(sz_plot, 700, 0, chooseSecond<RGB>);



    stbi::ImageWriter<stbi::format::PNG> writer;


    writer.writeInto("./image/spin/spin_evo", spin_img[0], spin_img.ncols(), spin_img.nrows());


	std::cout << "\n>..SpinTest: end" << std::endl;
	BENCHMARK_END(spin_test);
}