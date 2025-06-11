#pragma once
#include "sim.h"

#include <numer_complex.h>
#include <numer_mat.h>



using CGLESim_Profile = struct {
	unsigned			serial;
	unsigned			frame_interval;
	double				dt;
	unsigned			N1;
	unsigned			N2;
	double				L1;
	double				ep;
	double				b;
	double				c;
	double				d;
	double				e;
};


class CGLESim : public Sim
{
private:
	numer::mat<numer::Complex> field_;
	CGLESim_Profile profile_;


public:
	CGLESim(unsigned Steps);
	numer::mat<numer::RGB> renderFrame() override;
	~CGLESim() override;
};

