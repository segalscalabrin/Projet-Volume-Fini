#ifndef MAIN_LOOP_H
#define MAIN_LOOP_H

#include "Neos.hpp"
#include <vector>
#include <string>

using namespace neos;

void timeLoop(Grid *g, double final_time, int time_steps, double dt);

void computeLevelSet(Grid *g, double moving_line_y, std::vector<double> &levelset_values);
void applyRefinement(Grid *g, double moving_line_y);
void exportResults(Grid *g, int step, std::vector<double> &levelset_values);

#endif // MAIN_LOOP_H