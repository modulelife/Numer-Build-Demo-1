#pragma once
#include <numer_mat.h>
#include <numer_visualize.h>

class Sim
{
protected:
	unsigned step_cur_{ 0 };
	unsigned steps_;

public:
	Sim(unsigned Steps) : steps_(Steps) {}
	virtual ~Sim() = default;

	virtual numer::mat<numer::RGB> renderFrame() = 0;
};

