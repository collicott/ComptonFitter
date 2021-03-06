#include <iostream>
#include <APLCON.hpp>
#include "CEFTwrapper.h"
#include "DataHandling.h"
#include "TFile.h"
#include "TNtuple.h"

using namespace std;

struct params {
    double alpha;
    double beta;
    double E1E1;
    double M1M1;
    double E1M2;
    double M1E2;

    params(double _alpha,
           double _beta,
           double _E1E1,
           double _M1M1,
           double _E1M2,
           double _M1E2) :
        alpha(_alpha),
        beta(_beta),
        E1E1(_E1E1),
        M1M1(_M1M1),
        E1M2(_E1M2),
        M1E2(_M1E2)
    {}
};

struct constraints {
    double value;
    double error;

    constraints(double _value,
           double _error) :
        value(_value),
        error(_error)
    {}
};

int main(int argc, char *argv[])
{
    // Set up APLCON
    APLCON::Fit_Settings_t settings = APLCON::Fit_Settings_t::Default;
    settings.MaxIterations = 1000;
    settings.DebugLevel = 4;
    APLCON aplcon("Polarisabilities", settings);

	// Set up fitter
    Fitter myfit;

    // Set initial values
    params fitparam(10.2, 3.5, -4.0, 2.75, 0.22, 2.0);
    constraints alpha_beta_sum(13.8, 0.4);
    constraints alpha_beta_diff(7.6, 1.7);

    // Add params as unmeasured variables
    aplcon.AddUnmeasuredVariable("alpha",fitparam.alpha);
    aplcon.AddUnmeasuredVariable("beta", fitparam.beta);

    // Add some measured constraints
    aplcon.AddMeasuredVariable("alpha_beta_sum",alpha_beta_sum.value, alpha_beta_sum.error);
    aplcon.AddMeasuredVariable("alpha_beta_diff",alpha_beta_diff.value, alpha_beta_diff.error);

    // Set up constraint lamda functions
    auto sum_constraint  = [] (double a, double b, double sum)  {return sum - a - b;  };
    auto diff_constraint = [] (double a, double b, double diff) {return diff - a + b; };

    // Apply the constraints
//    aplcon.AddConstraint("Baldin", {"alpha", "beta", "alpha_beta_sum"},  sum_constraint);
    aplcon.AddConstraint("Diff",   {"alpha", "beta", "alpha_beta_diff"}, diff_constraint);

    // Set up Data handler and add data points as measured variables
    DataHandling dataHandler;
    dataHandler.process_file_list(argv[1]);
    vector<data> datapoint = dataHandler.GetDataList();

    for (auto i = 0; i < datapoint.size(); i++)
    {
        stringstream ss;
        ss << "datapoint" << i;
        aplcon.AddMeasuredVariable(ss.str(),
                                   datapoint[i].observable,
                                   datapoint[i].error);
    }

    // Create equality constraint for datapoints
    for (auto i = 0; i < datapoint.size(); i++)
    {
        stringstream constraintname;
        constraintname << "datapoint_" << i;

        stringstream pointname;
        pointname << "datapoint" << i;

        // setup a lambda function which returns 0
        auto equality_constraint = [&myfit,&datapoint,i,&fitparam]
                (double experiment,
                double alpha,
                double beta)
        {
            auto theory = myfit.Fit(datapoint[i].theta,
                                    datapoint[i].energy,
                                    alpha,
                                    beta,
                                    fitparam.E1E1,
                                    fitparam.M1M1,
                                    fitparam.E1M2,
                                    fitparam.M1E2);

            if (datapoint[i].data_type == Sigma_3)
                return experiment - theory.GetSigma3();

            else if (datapoint[i].data_type == Sigma_2x)
                return experiment - theory.GetSigma2x();

            else if (datapoint[i].data_type == Sigma_2z)
                return experiment - theory.GetSigma2z();

            else if (datapoint[i].data_type == Cross)
                return experiment - theory.GetCross();

            else cout << "ERROOOROROROROROROR" << endl;
            return 1000000.0;
        };

        aplcon.AddConstraint(
                    constraintname.str(),
                    {pointname.str(),
                     "alpha",
                     "beta"},
                     equality_constraint);
    }

    // do the fit, obtain ra structure
    const APLCON::Result_t& ra = aplcon.DoFit();
    cout << ra << endl;

    cout << "Calls to Executable: " << myfit.GetNCalls() << endl;

    return 0;
}
