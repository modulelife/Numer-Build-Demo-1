#pragma once
#include "sim.h"
#include <vector>
#include <numer_complex.h>
#include <numer_mat.h>
#include <numer_visualize.h>
#include <numer_qmkit.h>


class HO1dSim : public Sim
{
private:
	std::vector<std::vector<numer::Complex>> nbasis_bra_;
	std::vector<numer::Complex> nket_;
	numer::qm::HusimiQCalculator q_distri_gen_;

public:
	HO1dSim(unsigned Steps);
	numer::mat<numer::RGB> renderFrame() override;

};
