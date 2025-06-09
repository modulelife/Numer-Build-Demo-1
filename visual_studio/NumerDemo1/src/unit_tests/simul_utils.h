#pragma once
#include <string>
#include <sstream>
#include <iomanip>

const char* root_defaut = "./image/simul/";
const char* image_dir = "frames/";


class SimImagePath
{
private:
	std::string sim_name_;
	unsigned fnum_begin_;

	SimImagePath() = delete;
	SimImagePath(SimImagePath const&) = delete;
	SimImagePath(SimImagePath&&) = delete;
	SimImagePath& operator=(SimImagePath const&) = delete;

public:

	SimImagePath(const std::string& Sim_name, unsigned Frame_number_begin_)
		: sim_name_(Sim_name), fnum_begin_(Frame_number_begin_)
	{}


	const std::string getPath(unsigned count) const
	{
		std::ostringstream oss;
		oss << std::setfill('0') << std::setw(6) << count + fnum_begin_;
		return
			std::string(root_defaut) +
			sim_name_ + "/" +
			std::string(image_dir) +
			sim_name_ +
			oss.str();
	}
};