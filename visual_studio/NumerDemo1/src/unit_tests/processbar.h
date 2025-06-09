#pragma once
#include <iostream>
#include <sstream>
#include <atomic>
#include <chrono>
#include <thread>


class ProcessBar
{
private:
	enum { BAR_LENGTH = 44 };

	std::atomic<unsigned int>& _counter;
	unsigned int _total;

	ProcessBar() = delete;
	ProcessBar(ProcessBar const&) = delete;
	ProcessBar(ProcessBar&&) = delete;
	ProcessBar& operator=(ProcessBar const&) = delete;

public:

	ProcessBar(std::atomic<unsigned int>& counter, unsigned int total)
		: _counter(counter), _total(total)
	{}

	void track(const char* title)
	{
		std::cout << title << '\n';
		unsigned cnt = 0;
		while (cnt < _total)
		{
			std::ostringstream cmd_line;
			cnt = _counter.load(std::memory_order_relaxed);
			cmd_line << '\r';
			unsigned i = 0;
			for (; i < (cnt * BAR_LENGTH / _total); i++) cmd_line << '#';
			for (; i < BAR_LENGTH; i++) cmd_line << '.';
			cmd_line << "| " << cnt * 100 / _total << '%';
			cmd_line << " " << _total << " in total";
			std::cout << cmd_line.str();
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
		}
		std::cout << '\r';
		for (unsigned j = 0; j < BAR_LENGTH; j++) std::cout << '#';
		std::cout << "| done" << std::endl;
	}

	void track_period(const char* title, std::atomic<float>& sec_per_period, const char* item)
	{
		std::cout << title << '\n';
		unsigned cnt = 0;
		while (cnt < _total)
		{
			std::ostringstream cmd_line;
			cnt = _counter.load(std::memory_order_relaxed);
			cmd_line << '\r';
			unsigned i = 0;
			for (; i < (cnt * BAR_LENGTH / _total); i++) cmd_line << '#';
			for (; i < BAR_LENGTH; i++) cmd_line << '.';
			cmd_line << "| " << cnt * 100 / _total << '%';
			cmd_line << " " << sec_per_period.load(std::memory_order_relaxed) << item;
			cmd_line << " " << _total << " in total";
			std::cout << cmd_line.str();
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
		}
		std::cout << '\r';
		for (unsigned j = 0; j < BAR_LENGTH; j++) std::cout << '#';
		std::cout << "| done --.---" << item << " " << _total << " in total" << std::endl;
	}

};



static void monitor_task__(const char* title, std::atomic<unsigned int>& counter, std::atomic<float>& timer, unsigned int total)
{
	std::cout << std::endl;
	ProcessBar bar(counter, total);
	bar.track_period(title, timer, "sec/frame");
}

//shutoff any other console IO between 'Tracker tk();' and 'tk.stop();'.
//use .stop() before following console IO to prevent funny output.
class Tracker
{
private:
	std::jthread t_;

	Tracker() = delete;
	Tracker(Tracker const&) = delete;
	Tracker(Tracker&&) = delete;
	Tracker& operator=(Tracker const&) = delete;

public:
	Tracker(const char* title, std::atomic<unsigned int>& counter, std::atomic<float>& timer, unsigned int total)
		:t_(
			monitor_task__,
			title,
			std::ref(counter),
			std::ref(timer),
			total
		)
	{}

	void stop()
	{
		while (!t_.joinable());
		t_.join();
	}
};