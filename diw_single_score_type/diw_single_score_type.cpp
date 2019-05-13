// full_run_subarray_sum.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include <ctime>
#include <map>
#include <vector>
#include <ratio>
#include <windows.h>
#include <stdlib.h>
#include <fcntl.h>
#include <queue>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <functional>
#include <numeric>
#include <tuple>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/Peak1D.h>

#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>


#include <boost/tokenizer.hpp>

enum ScoreType { ALL_ADDED, RATIO, TARGET_INT };

class WindowSpec {

public:
	WindowSpec(int s, std::string w, std::string c) : sister(s), width(w), center(c) {};

	int sister;
	std::string width;
	std::string center;




	friend bool operator>(const WindowSpec &w1, const WindowSpec &w2);
	friend bool operator<(const WindowSpec &w1, const WindowSpec &w2);
	friend bool operator==(const WindowSpec &w1, const WindowSpec &w2);


};

class Offset {

public:
	Offset(int s, double o) : sister(s), offset(o) {};

	int sister;
	double offset;

};

class peak
{
public:
	double mz;
	double intensity;

	peak(double m = 0.0, double i = 0.0) : mz(m), intensity(i) {};
};

bool operator> (const WindowSpec &w1, const WindowSpec &w2)
{
	if (w1.center > w2.center)
		return true;
	else if (w1.center == w2.center)
		return w1.width > w2.width;
	else
		return false;
}

bool operator< (const WindowSpec &w1, const WindowSpec &w2)
{
	if (w1.center < w2.center)
		return true;
	else if (w1.center == w2.center)
		return w1.width < w2.width;
	else
		return false;
}

bool operator== (const WindowSpec &w1, const WindowSpec &w2)
{
	//bool center_close = (abs(atof(w1.center.c_str()) - atof(w2.center.c_str())) <= 0.01);
	//bool off_close = (abs(atof(w1.width.c_str()) - atof(w2.width.c_str())) <= 0.01);
	return w1.center == w2.center && w1.width == w2.width;
}

//
//	A lot of commented stuff in here is legacy code and/or an old version of the algorithm and should not be regarded seriously 
//
double windowScore(OpenMS::MSSpectrum<>& ms1, double mono_mz, int charge, double range_front, double range_end, int& skip_mono, int& lr, int &rr, std::ofstream& missing_f, int scan_num, ScoreType type)
{
	const double massDiff = 1.0033548378;

	std::vector<peak> peaks;
	std::vector<peak> light_spec;
	std::vector<peak> iso_hits;

	OpenMS::IsotopeDistribution id(5);
	id.estimateFromPeptideWeight(mono_mz * charge);
	id.renormalize();

	id.getContainer().at(0).first;
	
	for (auto it = ms1.begin(); it != ms1.end(); it++)
	{
		peaks.push_back(peak(it->getPos(), it->getIntensity()));
	}

	double target_isos[5] = { mono_mz, mono_mz + (massDiff / charge), mono_mz + 2 * (massDiff / charge), mono_mz + 3 * (massDiff / charge), mono_mz + 4 * (massDiff / charge) };
	double target_intensities[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

	int temp_iso = 0;
	while (temp_iso < 5 && target_isos[temp_iso] < range_front)
		++temp_iso;

	if (temp_iso == 5 || target_isos[temp_iso] > range_end)
	{
		std::cout << "Dyn window " << scan_num << " contains no isotope" << std::endl;
		std::cout << ms1.size() << " " << mono_mz << " " << charge << " " << range_front << " " << range_end << std::endl;
	}
	
	double tol = 20 * mono_mz * charge * (1 / 1000000.0);
	
	std::vector<peak>::iterator copier = peaks.begin();
	std::vector<peak>::iterator look_ahead = peaks.begin();

	int curr_iso = 0;
	bool hit_iso = false;
	double iso_comp = target_isos[0];

	while (look_ahead != peaks.end() && look_ahead->mz < range_front)
	{
		look_ahead++;
		copier++;
	}

	if (look_ahead == peaks.end())
		return -1;

	if (mono_mz < range_front && abs(target_isos[0] - copier->mz) >= tol)
	{
		++skip_mono;
		
		missing_f << scan_num << std::endl;
		
		copier = peaks.begin();
		while (copier != peaks.end() && abs(target_isos[0] - copier->mz) >= tol)
			copier++;

		double mono_diff = abs(target_isos[0] - copier->mz);
		double mono_last = copier->intensity;
		while (abs(target_isos[0] - copier->mz) < tol)
		{
			double mono_curr_diff = abs(target_isos[0] - copier->mz);
			if (mono_curr_diff < mono_diff)
			{
				mono_diff = mono_curr_diff;
				mono_last = copier->intensity;
			}
			++copier;
		}

		target_intensities[0] = mono_last;
		for (int other_iso = 1; other_iso < 5; ++other_iso)
		{
			target_intensities[other_iso] = target_intensities[0] * (id.getContainer().at(other_iso).second / id.getContainer().at(0).second);
		}
		curr_iso++;

		copier = peaks.begin();
		while (copier->mz < range_front)
			++copier;
	}

	double best_diff = 100000.0;
	double best_old_inten = 0.0;


	double total_pos = 0.0;
	double total_neg = 0.0;

	double measured_intensity[5] = { target_intensities[0], 0, 0, 0, 0 };

	while (copier != peaks.end() && copier->mz < target_isos[4] + 0.1 && copier->mz < range_end)
	{
		if (curr_iso < 5)
		{
			double curr_diff = abs(target_isos[curr_iso] - copier->mz);

			if (curr_diff < tol)
			{
				if (curr_diff < best_diff)
				{
					if (best_diff < 99999.0)
					{
						light_spec.at(light_spec.size() - 1) = peak(light_spec.at(light_spec.size() - 1).mz, -1 * best_old_inten);
						if (target_intensities[curr_iso] > best_old_inten)
						{
							total_pos -= best_old_inten;
							//std::cout << "Running Pos Count: " << total_pos << std::endl;
						}
						else
						{
							total_pos -= target_intensities[curr_iso];
							total_neg -= (best_old_inten - target_intensities[curr_iso]);
							//std::cout << "Running Pos Count: " << total_pos << std::endl;
							//std::cout << "Running Neg Count: " << total_neg << std::endl;
						}

						total_neg += best_old_inten;
						//std::cout << "Running Neg Count: " << total_neg << std::endl;

					}
					best_diff = curr_diff;
					best_old_inten = copier->intensity;


					if (curr_iso == 0)
						target_intensities[0] = best_old_inten;

					light_spec.push_back(peak(copier->mz, -1 * copier->intensity + 2 * target_intensities[curr_iso]));

					
					if (target_intensities[curr_iso] > copier->intensity)
					{
						total_pos += copier->intensity;
						//std::cout << "Running Pos Count: " << total_pos << std::endl;
					}
					else
					{
						total_pos += target_intensities[curr_iso];
						total_neg += (copier->intensity - target_intensities[curr_iso]);
						//std::cout << "Running Pos Count: " << total_pos << std::endl;
						//std::cout << "Running Neg Count: " << total_neg << std::endl;
					}

				}

				else
				{
					light_spec.push_back(peak(copier->mz, -1 * copier->intensity));
					total_neg += copier->intensity;
					//std::cout << "Running Neg Count: " << total_neg << std::endl;
				}
				copier++;

			}

			else if (copier->mz > target_isos[curr_iso])
			{
				if (curr_iso == 0)
				{
					for (int other_iso = 1; other_iso < 5; ++other_iso)
					{
						target_intensities[other_iso] = target_intensities[0] * (id.getContainer().at(other_iso).second / id.getContainer().at(0).second);
					}
				}

				++curr_iso;
				best_diff = 100000.0;
				best_old_inten = 0.0;


			}

			else
			{
				light_spec.push_back(peak(copier->mz, -1 * copier->intensity));
				total_neg += copier->intensity;
				//std::cout << "Running Neg Count: " << total_neg << std::endl;
				++copier;
			}

		}

		else if (copier->mz < range_end)
		{
			light_spec.push_back(peak(copier->mz, -1 * copier->intensity));
			total_neg += copier->intensity;
			//std::cout << "Running Neg Count: " << total_neg << std::endl;
			++copier;
		}

	}

	while (copier != peaks.end() && copier->mz < range_end)
	{
		light_spec.push_back(peak(copier->mz, -1 * copier->intensity));
		total_neg += copier->intensity;
		//std::cout << "Running Neg Count: " << total_neg << std::endl;
		copier++;
	}


	copier = peaks.begin();
	while (copier->mz < target_isos[0] - 0.4)
		++copier;

	--copier;
	curr_iso = 0;
	best_diff = 100000.0;
	best_old_inten = 0.0;

	while (copier != peaks.end() && copier->mz < target_isos[4] + 0.1)
	{
		if (curr_iso < 5)
		{
			double curr_diff = abs(target_isos[curr_iso] - copier->mz);

			//std::cout << "Target: " << target_isos[curr_iso] << std::endl;
			//std::cout << "Current: " << copier->mz << std::endl;
			//std::cout << "Diff: " << curr_diff << std::endl;

			if (curr_diff < tol)
			{
				if (curr_diff < best_diff)
				{
					best_diff = curr_diff;
					best_old_inten = copier->intensity;

					measured_intensity[curr_iso] = best_old_inten;

					if (curr_iso == 0)
						target_intensities[0] = best_old_inten;
					//if (target_intensities[curr_iso] > copier->intensity && abs((range_end - range_front) - 1.6) > .01)
					//	std::cout << "candidate for testing..." << std::endl

				}

				copier++;

			}

			else if (copier->mz > target_isos[curr_iso])
			{
				if (curr_iso == 0)
				{
					for (int other_iso = 1; other_iso < 5; ++other_iso)
					{
						target_intensities[other_iso] = target_intensities[0] * (id.getContainer().at(other_iso).second / id.getContainer().at(0).second);
					}
				}

				++curr_iso;
				best_diff = 100000.0;
				best_old_inten = 0.0;


			}

			else
			{
				//std::cout << "Running Neg Count: " << total_neg << std::endl;
				++copier;
			}

		}

		else
		{
			++copier;
		}

	}

	//standard total score
	double total_score = 0.0;

	for (int k = 0; k < light_spec.size(); k++)
	{
		total_score += light_spec.at(k).intensity;
	}

	if (type == ALL_ADDED)
		return total_score;

	else if (type == RATIO)
	{
		double ratio = total_pos / (total_neg + total_pos);

		//std::cout << "Ratio = " << ratio << std::endl;

		if (abs(total_neg + total_pos - 0.0) < 0.0001)
			return -1.0;

		return ratio;
	}

	else if (type == TARGET_INT)
		return total_pos;

	else
		return -1;

}


class targetScore
{
public:
	double target;
	double score;

	targetScore(double t = 0.0, double s = 0.0) : target(t), score(s) {};
};

class family
{
public:
	int parent;
	int sister;

	family(int p = 0, int s = 0) : parent(p), sister(s) {};
	friend bool operator<(const family &f1, const family &f2);
};

bool operator< (const family &f1, const family &f2)
{
	if (f1.parent < f2.parent)
		return true;
	else if (f1.parent == f2.parent)
		return f1.sister < f2.sister;
	else
		return false;
}

int main(int argc, char * argv[])
{

	std::ifstream dyn_file(argv[1]);

	std::ifstream window_file(argv[2]);

	std::ofstream log("log.txt");

	std::map<int, bool> dyn_scans;

	std::map<WindowSpec, std::vector<Offset>> offsets;

	if (dyn_file)
	{
		std::string line;
		int scan;

		while (std::getline(dyn_file, line))
		{
			std::istringstream iss(line);
			iss >> scan;

			dyn_scans.insert({ scan, true });

		}

	}

	if (window_file)
	{
		std::string line;
		int scan;
		std::string width;
		std::string center;
		double offset;

		while (std::getline(window_file, line))
		{
			std::istringstream iss(line);
			iss >> scan >> center >> offset >> width;

			center = center.substr(0, center.find('.') + 3);
			width = width.substr(0, width.find('.') + 2);

			log << "Adding " << scan << " " << width << " " << center << std::endl;


			auto offset_finder = offsets.find(WindowSpec(scan, width, center));
			if (offset_finder == offsets.end())
			{
				std::vector<Offset> temp;
				temp.push_back(Offset(scan, offset));
				offsets.insert({ WindowSpec(scan, width, center), temp });
			}

			else
			{
				offset_finder->second.push_back(Offset(scan, offset));
			}

		}

	}

	std::map<WindowSpec, std::string> dyn_parents;

	std::map<std::string, std::vector<family>> offset_family_lists;

	std::map<int, OpenMS::MSSpectrum<>> ms1_collection;

	std::map<int, double> ms2_scores;

	std::queue<targetScore> regular_scores;

	int rounding_prec = 5;
	int dynamic_less = 0;
	int dynamic_more = 0;
	int dyn_skip = 0;
	int reg_skip = 0;

	int left_rounding = 0;
	int right_rounding = 0;

	std::string score_type = argv[4];
	ScoreType scoring_type;

	if (score_type == "RATIO")
		scoring_type = RATIO;


	else if (score_type == "ALL_ADDED")
		scoring_type = ALL_ADDED;

	else
		scoring_type = ALL_ADDED;

	std::ofstream dyn_f;
	dyn_f.open(argv[5], std::ios_base::app);
	
	std::ofstream reg_f;
	reg_f.open(argv[6], std::ios_base::app);

	std::ofstream missing_f;
	missing_f.open("missing_mono.txt");

	std::ofstream missing_reg_f;
	missing_reg_f.open("missing_reg_mono.txt");

	OpenMS::MzMLFile mzMLDataFileProfile;
	OpenMS::MSExperiment msExperimentProfile;

	try {
		mzMLDataFileProfile.load(argv[3], msExperimentProfile);
	}
	catch (std::exception& e) {
		std::cout << e.what() << std::endl;
		//usage();
		return 1;
	}

	OpenMS::MSSpectrum<> s;
	int ms1_count = 0;
	int reg_ms2_count = 0;
	int dyn_ms2_count = 0;
	int found_dynamic = 0;
	int missing_dynamic = 0;

	int list_scan = 0;
	int diff = 0;

	std::string scan;

	int previous_ms1 = 1;

	for (int i = 0; i < msExperimentProfile.getNrSpectra(); i++)
	{
		s = msExperimentProfile.getSpectrum(i);

		std::string junk;
		std::istringstream natID = std::istringstream(s.getNativeID());
		natID >> junk >> junk >> scan;
		scan = scan.substr(scan.find('=') + 1);

		int scan_num = atoi(scan.c_str());

		if (s.getMSLevel() == 1)
		{
			std::string junk;
			std::istringstream natID = std::istringstream(s.getNativeID());
			natID >> junk >> junk >> scan;
			scan = scan.substr(scan.find('=') + 1);

			int scan_num = atoi(scan.c_str());

			previous_ms1 = scan_num;

			ms1_collection.insert({ scan_num, s });

		}

		else
		{
			double inject_time;
			
			auto scan_finder = dyn_scans.find(scan_num);
			if (scan_finder != dyn_scans.end())
			{
				dyn_ms2_count++;

				auto huh = s.getPrecursors().at(0).getPos();

				double window_width = 2 * s.getPrecursors().at(0).getIsolationWindowUpperOffset();
				std::string width_s = std::to_string(window_width);
				width_s = width_s.substr(0, width_s.find('.') + 2);

				std::string center_s = std::to_string(huh);
				center_s = center_s.substr(0, center_s.find('.') + 3);

				bool look_lower = false;
				bool look_higher = false;
				bool found = false;

				double offset;
				int sister;

				auto window_finder = offsets.find(WindowSpec(scan_num, width_s, center_s));
				if (window_finder != offsets.end())
				{
					//std::cout << "Things as expected" << std::endl;
					int distance = 9999;
					for (int i = 0; i < window_finder->second.size() && window_finder->second.at(i).sister < scan_num; i++)
					{
						distance = scan_num - window_finder->second.at(i).sister;
						offset = window_finder->second.at(i).offset;

						sister = window_finder->second.at(i).sister;
					}
					if (distance > 50)
						look_lower = true;

					else
						found = true;
				}

				else
				{
					look_lower = true;
				}

				if (look_lower)
				{
					double temp = huh - .01;
					center_s = std::to_string(temp);
					center_s = center_s.substr(0, center_s.find('.') + 3);

					auto window_finder = offsets.find(WindowSpec(scan_num, width_s, center_s));
					if (window_finder != offsets.end())
					{
						std::cout << "Things as expected" << std::endl;
						int distance = 9999;
						for (int i = 0; i < window_finder->second.size() && window_finder->second.at(i).sister < scan_num; i++)
						{
							distance = scan_num - window_finder->second.at(i).sister;
							offset = window_finder->second.at(i).offset;

							sister = window_finder->second.at(i).sister;

						}
						if (distance > 50)
							look_higher = true;

						else
							found = true;
					}

					else
					{
						look_higher = true;
					}

				}

				if (look_higher)
				{
					double temp = huh + .01;
					center_s = std::to_string(temp);
					center_s = center_s.substr(0, center_s.find('.') + 3);

					auto window_finder = offsets.find(WindowSpec(scan_num, width_s, center_s));
					if (window_finder != offsets.end())
					{
						std::cout << "Things as expected" << std::endl;
						int distance = 9999;
						for (int i = 0; i < window_finder->second.size() && window_finder->second.at(i).sister < scan_num; i++)
						{
							distance = scan_num - window_finder->second.at(i).sister;
							offset = window_finder->second.at(i).offset;

							sister = window_finder->second.at(i).sister;
						}
						if (distance > 50)
							look_higher = true;

						else
							found = true;
					}

					else
					{
						look_higher = true;
					}
				}

				if (!found)
				{
					log << "Issues with scan: " << scan_num << " " << width_s << " " << center_s << std::endl;
				}

				double scan_offset = offset;
				double final_center = s.getPrecursors().at(0).getPos() - scan_offset;


				double mono_mz = (double)(s.getPrecursors().at(0).getMZ() - scan_offset);
				int charge = s.getPrecursors().at(0).getCharge();
				double front = (double)final_center - s.getPrecursors().at(0).getIsolationWindowLowerOffset();
				double end = (double)final_center + s.getPrecursors().at(0).getIsolationWindowUpperOffset();

				//this is to get data about edge conditions
				double score= 0.0;

				if (ms1_collection.find(previous_ms1) != ms1_collection.end())
				{
					score = windowScore(ms1_collection.find(previous_ms1)->second, mono_mz, charge, front, end, dyn_skip, left_rounding, right_rounding, missing_f, scan_num, scoring_type);	
				}

				dyn_f << score << std::endl;
				
				auto sister_score_finder = ms2_scores.find(sister);
				if (sister_score_finder != ms2_scores.end())
				{
					if (score < sister_score_finder->second)
					{
						dynamic_less++;
					}

					else
						dynamic_more++;
				}

			}

			else
			{
				std::vector<OpenMS::String> keyss;

				s.getKeys(keyss);

				std::string scan_spec;

				try
				{
					scan_spec = s.getMetaValue(keyss.at(6));
				}

				catch (...)
				{
					std::cout << "Issue with keys, breaking out of scan" << std::endl;
					break;
				}

				std::vector<std::string> scan_tokens;

				boost::char_separator<char> sep(" ");
				boost::tokenizer<boost::char_separator<char>> tokens(scan_spec, sep);
				for (const auto& t : tokens) {
					scan_tokens.push_back(t);
				}

				double final_center = 0.0;

				for (int k = 0; k < scan_tokens.size(); k++)
				{
					std::string curr_tok = scan_tokens.at(k);

					if (curr_tok.find('@') != std::string::npos)
					{
						std::string center_str = curr_tok.substr(0, curr_tok.find('@'));
						final_center = atof(center_str.c_str());
						break;
					}
				}


				double mono_mz = (double)s.getPrecursors().at(0).getMZ();
				int charge = s.getPrecursors().at(0).getCharge();
				double front = (double)final_center - s.getPrecursors().at(0).getIsolationWindowLowerOffset();
				double end = (double)final_center + s.getPrecursors().at(0).getIsolationWindowUpperOffset();

				double score = 0.0;

				if (ms1_collection.find(previous_ms1) != ms1_collection.end())
				{
					score = windowScore(ms1_collection.find(previous_ms1)->second, mono_mz, charge, front, end, reg_skip, left_rounding, right_rounding, missing_reg_f, scan_num, scoring_type);
				}
				
				reg_f << score << std::endl;

				regular_scores.push(targetScore(mono_mz, score));
				ms2_scores.insert({ scan_num, score });

				reg_ms2_count++;
			}

		}


		if (i % 100 == 0)
		{
			std::cout << i << " scans done" << std::endl;

			std::map<int, OpenMS::MSSpectrum<>>::iterator ms1_it;
			for (ms1_it = ms1_collection.begin(); ms1_it != ms1_collection.end(); ms1_it++)
			{
				if ((i - ms1_it->first) > 60)
				{
					ms1_collection.erase(ms1_it);
					ms1_it = ms1_collection.begin();
				}
			}
		}

	}

	//close files and report numbers
	dyn_f.close();
	reg_f.close();
	missing_f.close();
	missing_reg_f.close();

	std::ifstream scores("temp_scores.txt");

	std::cout << "Found Dyn = " << found_dynamic << " Missing Dyn = " << missing_dynamic << std::endl;
	std::cout << dyn_ms2_count << std::endl;
	std::cout << reg_ms2_count << std::endl;
	std::cout << "Dynamic More: " << dynamic_more << " Dynamic Less: " << dynamic_less << std::endl;
	std::cout << dyn_skip << " Dyn windows skip monoisotope" << std::endl;
	std::cout << reg_skip << " Reg windows skip monoisotope" << std::endl;
	std::cout << "Left rounding issues: " << left_rounding << std::endl;
	std::cout << "Right rounding issues: " << right_rounding << std::endl;
	std::cout << "Finito" << std::endl;

	return 0;
}

