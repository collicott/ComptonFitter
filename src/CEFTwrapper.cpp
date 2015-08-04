#include "CEFTwrapper.h"

#include <iostream>
#include <stdexcept>

using namespace std;

double CEFTwrapper::read_double(const string &filename) const
{
    const string absolute_filename(folder+"/"+filename);

    ifstream input;

    input.open(absolute_filename);

    if(!input.is_open())
        throw result_file_error(absolute_filename);

    double dummy=0.0;
    double value=0.0;
    input >> dummy >> value;
    input.close();
    return value;
}
//tter::~CEFTwrapper() {
//

double CEFTwrapper::Fit(const double th, const double E, const double a, const double b, const double E1E1, const double M1M1, const double E1M2, const double M1E2, std::string data_type)
{
    ++n_calls;

    string params =  to_string(th) + "_" + to_string(E) + "_" + to_string(a) + "_" + to_string(b) + "_" + to_string(E1E1) + "_" + to_string(M1M1) + "_" + to_string(E1M2) + "_" + to_string(M1E2);
    string map = data_type + "_" + params;

    if(results_map.find(map) != results_map.end())
    {
        // point already calculated, no need to recalculate
        n_folders++;
    }
    else
    {
        // Create unique folder
        char tmp[256] = "/tmp/comptonfit_XXXXXX";
        folder = mkdtemp(tmp);

        // Throw error if
        if(folder.empty())
            throw std::runtime_error("Can't create CEFTwrapper output directory");

        // Create command, pass parameters and folder name
        string cmd = command + " "
                + to_string(th) + " " + to_string(E) + " " + to_string(a) + " " + to_string(b) + " " + to_string(E1E1) + " " + to_string(M1M1) + " " + to_string(E1M2) + " " + to_string(M1E2) + " " + folder;

        // Run command
        int r = system(cmd.c_str());

        // Throw error if command fails
        if(r!=0) throw fit_command_error(r);

        // Make a new map entry for data points
        string sigma2x = "Sigma_2x_"+ params;
        results_map[sigma2x] = GetSigma2x();

        string sigma2z = "Sigma_2z_"+ params;
        results_map[sigma2z] = GetSigma2z();

        string sigma3  = "Sigma_3_"+ params;
        results_map[sigma3] = GetSigma3();

        string cross   = "Cross_"+ params;
        results_map[cross] = GetCross();

        // Remove folder
        system(string("rm -rf "+folder).c_str());
    }

    return results_map[map];
}

double CEFTwrapper::GetSigma2x() const
{
    return read_double("xSigma2x_Adistribution.output");
}

double CEFTwrapper::GetSigma2z() const
{
    return read_double("xSigma2z_Adistribution.output");
}

double CEFTwrapper::GetSigma3() const
{
    return (read_double("xBa_Adistribution.output"))/100.0;
}

double CEFTwrapper::GetCross() const
{
    return read_double("xcs_Adistribution.output");
}
