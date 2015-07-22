#include "CEFTwrapper.h"

#include <iostream>
#include <stdexcept>

using namespace std;

double Fitter::Result::read_double(const string &filename) const
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

Fitter::Result::~Result() {
  //  if(!folder.empty()) {
//        cout << "removing " + folder << endl;
   //     system(string("rm -rf "+folder).c_str());
  //  }
}

Fitter::Result Fitter::Fit(const double th, const double E, const double a, const double b, const double E1E1, const double M1M1, const double E1M2, const double M1E2) {

    ++n_calls;

    // Make folder
    string folder = "tmp/Pascalutsa_" + to_string(th) + "_" + to_string(E) + "_" + to_string(a) + "_" + to_string(b) + "_" + to_string(E1E1) + "_" + to_string(M1M1) + "_" + to_string(E1M2) + "_" + to_string(M1E2);

    // Check if folder exists already
    struct stat sb;
    if (stat(folder.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
    {
        n_folders++;
       // folder exists, return result without recalculating
       return std::move(Result(folder));
    }
    // else continue, and create folder, run command.. jadda jadda
    else
        system(string("mkdir "+folder).c_str());

    // Create command, pass parameters and folder name

    string cmd = command + " "
            + to_string(th) + " " + to_string(E) + " " + to_string(a) + " " + to_string(b) + " " + to_string(E1E1) + " " + to_string(M1M1) + " " + to_string(E1M2) + " " + to_string(M1E2) + " " + folder;

    // Run command
    int r = system(cmd.c_str());

    // Throw error
    if(r!=0) throw fit_command_error(r);

    return std::move(Result(folder));
}

double Fitter::Result::GetSigma2x() const
{
    return read_double("xSigma2x_Adistribution.output");
}

double Fitter::Result::GetSigma2z() const
{
    return read_double("xSigma2z_Adistribution.output");
}

double Fitter::Result::GetSigma3() const
{
    return (read_double("xBa_Adistribution.output"))/100.0;
}

double Fitter::Result::GetCross() const
{
    return read_double("xcs_Adistribution.output");
}
