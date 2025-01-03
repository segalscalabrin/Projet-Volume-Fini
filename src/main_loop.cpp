#include "main_loop.h"

void timeLoop(Grid *g, double final_time, int time_steps, double dt) {
    for (int step = 0; step <= time_steps; ++step) {
        double current_time = step * dt;
        double moving_line_y = 1.0 - 2.0 * current_time;

        applyRefinement(g, moving_line_y);

        std::vector<double> levelset_values(g->nbCells());
        computeLevelSet(g, moving_line_y, levelset_values);

        exportResults(g, step, levelset_values);
    }
}

void computeLevelSet(Grid *g, double moving_line_y, std::vector<double> &levelset_values) {
    for (auto &cell : g->getCells()) {
        const long &id = cell.getId();
        std::array<double, 3> ic;

        if (g->getCell(id).isInterior()) {
            ic = g->evalCellCentroid(id);
            double y = ic[1];
            levelset_values[id] = y - moving_line_y;
        }
    }
}

void applyRefinement(Grid *g, double moving_line_y) {
    for (auto &cell : g->getCells()) {
        const long &id = cell.getId();
        std::array<double, 3> ic;

        if (g->getCell(id).isInterior()) {
            ic = g->evalCellCentroid(id);
            double y = ic[1];
            double levelset = y - moving_line_y;

            if (fabs(levelset) < 0.1) {
                g->markCellForRefinement(id);
            } else if (fabs(levelset) < 0.5) {
                g->markCellForRefinement(id);
            } else if (fabs(levelset) < 1.0) {
                g->markCellForCoarsening(id);
            } else {
                g->markCellForCoarsening(id);
            }
        }
    }
    g->update(true, true);
}

void exportResults(Grid *g, int step, std::vector<double> &levelset_values) {
    std::string export_name = "solutions/levelset_t" + std::to_string(step);
    g->setExportName(export_name);
    g->addData("LevelSet", levelset_values);
    g->write();
}