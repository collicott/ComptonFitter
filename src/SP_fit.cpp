#include <iostream>
#include <APLCON.hpp>
#include "CEFTwrapper.h"
#include "DataHandling.h"

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
    settings.DebugLevel = 5;
    APLCON aplcon("Polarisabilities", settings);

    APLCON::Variable_Settings_t var_settings = APLCON::Variable_Settings_t::Default;

	// Set up fitter
    Fitter myfit;

    // Set initial values
    params fitparam(12.0, 1.9, -3.0, 3.5, -0.8, 2.0);
    constraints alpha_beta_sum(13.8, 1.8);
    constraints alpha_beta_diff(10, 1.8);
    constraints gamma_0(1.01, 0.08);
    constraints gamma_pi(8.0, 1.8);

    // Add params as unmeasured variables
    aplcon.AddUnmeasuredVariable("alpha",fitparam.alpha);
    aplcon.AddUnmeasuredVariable("beta", fitparam.beta);
    aplcon.AddUnmeasuredVariable("E1E1",fitparam.E1E1,var_settings);
    aplcon.AddUnmeasuredVariable("M1M1",fitparam.M1M1,var_settings);
    aplcon.AddUnmeasuredVariable("E1M2",fitparam.E1M2,var_settings);
    aplcon.AddUnmeasuredVariable("M1E2",fitparam.M1E2,var_settings);

    aplcon.AddMeasuredVariable("alpha_beta_sum",alpha_beta_sum.value, alpha_beta_sum.error);
    aplcon.AddMeasuredVariable("alpha_beta_diff",alpha_beta_diff.value, alpha_beta_diff.error);
    aplcon.AddMeasuredVariable("gamma_0",gamma_0.value, gamma_0.error);
    aplcon.AddMeasuredVariable("gamma_pi",gamma_pi.value, gamma_pi.error);

    // Set up constraint lamda functions
    auto sum_constraint  = [] (double a, double b, double sum)  {return sum - a - b;  };
    auto diff_constraint = [] (double a, double b, double diff) {return diff - a + b; };
    auto g0_constraint   = [] (double E1E1, double M1M1, double E1M2, double M1E2, double g0 )
    {
        return g0 + E1E1 + M1M1 + E1M2 + M1E2;
    };
    auto gpi_constraint  = [] (double E1E1, double M1M1, double E1M2, double M1E2, double gpi )
    {
        return gpi + E1E1 - M1M1 + E1M2 - M1E2;
    };

    // Apply the constraints
    aplcon.AddConstraint("Baldin", {"alpha", "beta", "alpha_beta_sum"},  sum_constraint);
    aplcon.AddConstraint("Diff",   {"alpha", "beta", "alpha_beta_diff"}, diff_constraint);
    aplcon.AddConstraint("gamma0", {"E1E1", "M1M1", "E1M2", "M1E2", "gamma_0"},  g0_constraint);
    aplcon.AddConstraint("gammapi",{"E1E1", "M1M1", "E1M2", "M1E2", "gamma_pi"}, gpi_constraint);


    // Set up Data handler and add data points as measured variables
    DataHandling dataHandler;
    string filelist(argv[1]);
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
        auto equality_constraint = [&myfit,&datapoint,i]
                (double experiment,
                double alpha,
                double beta,
                double E1E1,
                double M1M1,
                double E1M2,
                double M1E2)
        {
            auto theory = myfit.Fit(datapoint[i].theta,
                                    datapoint[i].energy,
                                    alpha,
                                    beta,
                                    E1E1,
                                    M1M1,
                                    E1M2,
                                    M1E2);
//            cout << alpha<< " "
//               << beta<<" "
//               << E1E1<<" "
//               << M1M1 << " "
//               << E1M2 << " "
//               << M1E2 << endl;

            return experiment - theory.GetSigma3();
        };

        aplcon.AddConstraint(
                    constraintname.str(),
                    {pointname.str(),
                     "alpha",
                     "beta",
                     "E1E1",
                     "M1M1",
                     "E1M2",
                     "M1E2"},
                     equality_constraint);
    }

    // setup a lambda function which returns 0
    auto alpha_constraint = [] (double alpha, double alpha_theory)
    {
        return alpha - alpha_theory;
    };


    // do the fit, obtain ra structure
    const APLCON::Result_t& ra = aplcon.DoFit();
    cout << ra << endl;

    cout << "Calls to Executable: " << myfit.GetNCalls() << endl;

    return 0;
}

/*
int main()
{
   // Set up fitter
   Fitter myfit;
   params fitparam(12.0, 1.9, -4.5, 3.0, -0.8, 2.0);

   for (double E1E1 = -5.5; E1E1 < 5.5; E1E1+=0.25)
   {
       auto theory = myfit.Fit(90,
                               287,
                               fitparam.alpha,
                               fitparam.beta,
                               E1E1,
                               fitparam.M1M1,
                               fitparam.E1M2,
                               fitparam.M1E2);

       cout << E1E1 << " " << theory.GetSigma2x() << endl;

   }

   return 0;

}
*/
