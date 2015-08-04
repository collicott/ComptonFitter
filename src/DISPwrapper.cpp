#include "DISPwrapper.h"

#include <iostream>
#include <stdexcept>

using namespace std;

double DISPwrapper::read_double(const string &filename) const
{
    const string absolute_filename(filename);

    ifstream input;

    input.open(absolute_filename);

    if(!input.is_open())
        throw result_file_error(absolute_filename);

    std::string dummy;
    double value=0.0;
    input >> dummy >> value;
    input.close();
    return value;
}
//tter::~DISPwrapper() {
//

double DISPwrapper::Fit(const double th, const double E, const double a, const double b, const double E1E1, const double M1M1, const double E1M2, const double M1E2, std::string data_type)
{
    ++n_calls;

    string params = to_string(th) + "_" + to_string(E) + "_" + to_string(a) + "_" + to_string(b) + "_" + to_string(E1E1) + "_" + to_string(M1M1) + "_" + to_string(E1M2) + "_" + to_string(M1E2);
    string map = data_type + "_" + params;

    if(results_map.find(map) != results_map.end())
    {
        // point already calculated, no need to recalculate
        n_folders++;
    }
    else
    {
        // Remove old file because pasquini appends
        system(string("rm pasquini.txt").c_str());

        // Create command, pass parameters
        string type_string;              // Figure out data type string for command
        if(data_type == "Sigma_3")       type_string = "--sigma3";
        else if(data_type == "Sigma_2x") type_string = "--sigma2x";
        else if(data_type == "Sigma_2z") type_string = "--sigma2z";
        else if(data_type == "Cross")    type_string = "--dcs-cm";

        string space = " ";

        string cmd = command + space + type_string + space
                + "--energy-lab=" + to_string(E) + space
                + "--theta-lab=" + to_string(th) + space
                + "--alpha="+ to_string(a) + space
                + "--beta=" + to_string(b) + space
                + "--ge1e1="+ to_string(E1E1) + space
                + "--gm1m1="+ to_string(M1M1) + space
                + "--ge1m2="+ to_string(E1M2) + space
                + "--gm1e2="+ to_string(M1E2) + space;

        // Run command
        int r = system(cmd.c_str());

        // Throw error if command fails
        if(r!=0) throw fit_command_error(r);

        // Make a new map entry for data points, flip sigma3 values
        double theory = read_double("pasquini.txt");
        if (data_type == "Sigma_3") results_map[map] = (-1.0)*theory;
        else results_map[map] = theory;
    }

    return results_map[map];
}

