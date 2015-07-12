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
    if(!folder.empty()) {
//        cout << "removing " + folder << endl;
        system(string("rm -rf "+folder).c_str());
    }
}

Fitter::Result Fitter::Fit(const double th, const double E, const double a, const double b, const double E1E1, const double M1M1, const double E1M2, const double M1E2) {

    ++n_calls;

    char folder[256] = "/tmp/comptonfit_XXXXXX";

    string tmp = mkdtemp(folder);

    if(tmp.empty())
        throw std::runtime_error("Can't create fitter output directory");

    string cmd = command + " "
            + to_string(th) + " " + to_string(E) + " " + to_string(a) + " " + to_string(b) + " " + to_string(E1E1) + " " + to_string(M1M1) + " " + to_string(E1M2) + " " + to_string(M1E2) + " " + tmp;

  //  cout << "Running " << cmd << endl;
    int r = system(cmd.c_str());

    if(r!=0)
        throw fit_command_error(r);

    return std::move(Result(tmp));
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
