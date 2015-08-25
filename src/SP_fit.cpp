#include <iostream>
#include <APLCON.hpp>
#include "CEFTwrapper.h"
#include "DISPwrapper.h"
#include "DataHandling.h"
#include "DisplayResults.h"
#include "TFile.h"
#include "TNtuple.h"
#include <ctime>
#include <string>
#include <stdlib.h>
#include <cmath>

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


class SP_fit {
public:
    SP_fit() {}
    APLCON::Result_t run(int argc, char *argv[], vector<data> datapoint, string theory_code, params fitparam, int instance, TFile *f, TNtuple *ntuple);
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
    if (argc == 1)
    {
        cout << "Example: SP_fit data/SP_fit_Collicott_Martel.dat --Pascalutsa out.root 5" << endl;
        return 1;
    }

    SP_fit fitter;
    params fitparam(11.2, 2.5, 0, 0, 0, 0);
    double scaler_spread = 2;
    double spin_spread = 6;

    // Set up Data handler
    DataHandling dataHandler;
    dataHandler.process_file_list(argv[1]);
    vector<data> datapoint = dataHandler.GetDataList();

    // Quickly parse through arguments to find theory and filelist
    string theory_code;

    string flag(argv[2]);
    if (flag == "--Pascalutsa") theory_code = "Pascalutsa";
    else if (flag == "--Pasquini") theory_code = "Pasquini";
    else {
        cout << "Error: unknown theory " << argv[2] << endl;
        cout << "Using default Pascalutsa" << endl;
        theory_code = "Pascalutsa";
    }

    // Set up some ROOT storage
    string name;
    if (argc > 3) name = argv[3];
    else name = "SP_fit.root";
    TFile *f = new TFile(name.c_str(),"RECREATE");
    TNtuple *ntuple = new TNtuple("SP_fit",
                                  "SP_fit",
                                  "instance:call:alpha:beta:E1E1:M1M1:E1M2:M1E2");

    TNtuple *start_values = new TNtuple("start_values",
                                  "start_values",
                                  "instance:alpha_start:beta_start:E1E1_start:M1M1_start:E1M2_start:M1E2_start");


    TNtuple *result = new TNtuple("result",
                                  "result",
                                  "instance:alpha:alpha_error:beta:beta_error:E1E1:E1E1_error:M1M1:M1M1_error:E1M2:E1M2_error:M1E2:M1E2_error:chi");

    // Get number of fits
    int n_fits = 1;
    if(argv[4]) n_fits = atoi(argv[4]);

    //Type of random number distribution
    uniform_real_distribution<double> alpha(fitparam.alpha - scaler_spread, fitparam.alpha + scaler_spread);  //(min, max)
    uniform_real_distribution<double> beta(fitparam.beta - scaler_spread, fitparam.beta + scaler_spread);  //(min, max)
    uniform_real_distribution<double> E1E1(fitparam.E1E1 - spin_spread, fitparam.E1E1 + spin_spread);  //(min, max)
    uniform_real_distribution<double> M1M1(fitparam.M1M1 - spin_spread, fitparam.M1M1 + spin_spread);  //(min, max)
    uniform_real_distribution<double> E1M2(fitparam.E1M2 - spin_spread, fitparam.E1M2 + spin_spread);  //(min, max)
    uniform_real_distribution<double> M1E2(fitparam.M1E2 - spin_spread, fitparam.M1E2 + spin_spread);  //(min, max)

    mt19937 rng_alpha; rng_alpha.seed(random_device{}());
    mt19937 rng_beta; rng_beta.seed(random_device{}());
    mt19937 rng_E1E1; rng_E1E1.seed(random_device{}());
    mt19937 rng_M1M1; rng_M1M1.seed(random_device{}());
    mt19937 rng_E1M2; rng_E1M2.seed(random_device{}());
    mt19937 rng_M1E2; rng_M1E2.seed(random_device{}());

    // generate N random fits.
    double alpha_result = 0.0;
    double beta_result = 0.0;
    double E1E1_result = 0.0;
    double M1M1_result = 0.0;
    double E1M2_result = 0.0;
    double M1E2_result = 0.0;

    double dalpha_result = 0.0;
    double dbeta_result = 0.0;
    double dE1E1_result = 0.0;
    double dM1M1_result = 0.0;
    double dE1M2_result = 0.0;
    double dM1E2_result = 0.0;
    for (int i=0; i<n_fits; i++)
    {
        params random_fitparam(alpha(rng_alpha), beta(rng_beta), E1E1(rng_E1E1), M1M1(rng_M1M1), E1M2(rng_E1M2), M1E2(rng_M1E2));
        APLCON::Result_t ra = fitter.run(argc,argv,datapoint,theory_code,random_fitparam,i,f,ntuple);

        start_values->Fill(i,
                     ra.Variables.at("alpha").Value.Before,
                     ra.Variables.at("beta").Value.Before,
                     ra.Variables.at("E1E1").Value.Before,
                     ra.Variables.at("M1M1").Value.Before,
                     ra.Variables.at("E1M2").Value.Before,
                     ra.Variables.at("M1E2").Value.Before
                     );

        result->Fill(i,
                     ra.Variables.at("alpha").Value.After,
                     ra.Variables.at("alpha").Sigma.After,
                     ra.Variables.at("beta").Value.After,
                     ra.Variables.at("beta").Sigma.After,
                     ra.Variables.at("E1E1").Value.After,
                     ra.Variables.at("E1E1").Sigma.After,
                     ra.Variables.at("M1M1").Value.After,
                     ra.Variables.at("M1M1").Sigma.After,
                     ra.Variables.at("E1M2").Value.After,
                     ra.Variables.at("E1M2").Sigma.After,
                     ra.Variables.at("M1E2").Value.After,
                     ra.Variables.at("M1E2").Sigma.After,
                     ra.ChiSquare
                     );

        alpha_result+=ra.Variables.at("alpha").Value.After;
        beta_result+=ra.Variables.at("beta").Value.After;
        E1E1_result+=ra.Variables.at("E1E1").Value.After;
        M1M1_result+=ra.Variables.at("M1M1").Value.After;
        E1M2_result+=ra.Variables.at("E1M2").Value.After;
        M1E2_result+=ra.Variables.at("M1E2").Value.After;

        dalpha_result+=ra.Variables.at("alpha").Sigma.After;
        dbeta_result+=ra.Variables.at("beta").Sigma.After;
        dE1E1_result+=ra.Variables.at("E1E1").Sigma.After;
        dM1M1_result+=ra.Variables.at("M1M1").Sigma.After;
        dE1M2_result+=ra.Variables.at("E1M2").Sigma.After;
        dM1E2_result+=ra.Variables.at("M1E2").Sigma.After;
    }

    // Compute average parameters
    alpha_result = alpha_result/double(n_fits);
    beta_result = beta_result/double(n_fits);
    E1E1_result = E1E1_result/double(n_fits);
    M1M1_result = M1M1_result/double(n_fits);
    E1M2_result = E1M2_result/double(n_fits);
    M1E2_result = M1E2_result/double(n_fits);

    dalpha_result = dalpha_result/double(n_fits);
    dbeta_result = dbeta_result/double(n_fits);
    dE1E1_result = dE1E1_result/double(n_fits);
    dM1M1_result = dM1M1_result/double(n_fits);
    dE1M2_result = dE1M2_result/double(n_fits);
    dM1E2_result = dM1E2_result/double(n_fits);

    // Display results
//    DisplayResults display;
//    display.ShowFitResults(f,datapoint,alpha_result,beta_result,E1E1_result,M1M1_result,E1M2_result,M1E2_result,theory_code);

    f->Write();

    // Output final, averaged results
    cout << "****************************************" << endl;
    cout << "Averaged results:" << endl;
    cout << "alpha: " << alpha_result << " +- " << dalpha_result << endl;
    cout << "beta: " << beta_result << " +- " << dbeta_result << endl;
    cout << "E1E1: " << E1E1_result << " +- " << dE1E1_result << endl;
    cout << "M1M1: " << M1M1_result << " +- " << dM1M1_result << endl;
    cout << "E1M2: " << E1M2_result << " +- " << dE1M2_result << endl;
    cout << "M1E2: " << M1E2_result << " +- " << dM1E2_result << endl;

    return 0;

}

APLCON::Result_t SP_fit::run(int argc, char *argv[], vector<data> datapoint, string theory_code, params fitparam, int instance, TFile *f, TNtuple *ntuple)
{
    // Set up some time
    clock_t begin = clock();

    // Set up fitters
    CEFTwrapper CEFTfit;
    DISPwrapper DISPfit;

    // Set up APLCON
    APLCON::Fit_Settings_t settings = APLCON::Fit_Settings_t::Default;
    settings.MaxIterations = 1000;
    settings.DebugLevel = 4;
    APLCON aplcon("Polarisabilities", settings);

    // Set initial values
//  params fitparam(12.1, 1.6, -3.3, 3.0, 0.2, 1.1);     // (nominal BxPT) -3.7
//  params fitparam(12.1, 1.6, -4.3, 2.9, -0.02, 2.2);   // (nominal HDPV) -3.7
    
    // new alpha-beta
//    params fitparam(11.2, 2.5, -3.3, 3.0, 0.2, 1.1);     // (nominal BxPT) -3.7

//    constraints alpha_beta_diff(10.5, 1.6);  // Old numbers
    constraints alpha_beta_diff(7.6, 1.7); // Griesshammer
    constraints alpha_beta_sum(13.8, 0.4);
    constraints gamma_0(-1.01, 0.08);
    constraints gamma_pi(8.0, 1.8);

    // Add params as unmeasured variables
    aplcon.AddUnmeasuredVariable("alpha",fitparam.alpha);
    aplcon.AddUnmeasuredVariable("beta", fitparam.beta);
    aplcon.AddUnmeasuredVariable("E1E1",fitparam.E1E1);
    aplcon.AddUnmeasuredVariable("M1M1",fitparam.M1M1);
    aplcon.AddUnmeasuredVariable("E1M2",fitparam.E1M2);
    aplcon.AddUnmeasuredVariable("M1E2",fitparam.M1E2);

    // Add some measured constraints
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

    // Add datapoints as Measured
    for (auto i = 0; i < datapoint.size(); i++)
    {
        stringstream ss;
        ss << "datapoint" << i;
        aplcon.AddMeasuredVariable(ss.str(),
                                   datapoint[i].observable,
                                   datapoint[i].error);
    }

    int call = 0;
    // Create equality constraint for datapoints
    for (auto i = 0; i < datapoint.size(); i++)
    {
        stringstream constraintname;
        constraintname << "datapoint_" << i;

        stringstream pointname;
        pointname << "datapoint" << i;

        // setup a lambda function which returns 0
        auto equality_constraint = [&CEFTfit,&DISPfit,&datapoint,i,&ntuple,&call,&instance,&theory_code]
                (double experiment,
                double alpha,
                double beta,
                double E1E1,
                double M1M1,
                double E1M2,
                double M1E2)
        {
            double theory = 0;
            if(theory_code == "Pascalutsa")
                 theory = CEFTfit.Fit(datapoint[i].theta,
                                      datapoint[i].energy,
                                      alpha,
                                      beta,
                                      E1E1,
                                      M1M1,
                                      E1M2,
                                      M1E2,
                                      datapoint[i].data_type);
            else if (theory_code == "Pasquini")
                 theory = DISPfit.Fit(datapoint[i].theta,
                                      datapoint[i].energy,
                                      alpha,
                                      beta,
                                      E1E1,
                                      M1M1,
                                      E1M2,
                                      M1E2,
                                      datapoint[i].data_type);

            ntuple->Fill(instance,
                         call,
                         alpha,
                         beta,
                         E1E1,
                         M1M1,
                         E1M2,
                         M1E2);
            call++;

    //        if (theory_code == "Pasquini")
  //          cout << experiment << " " << theory << endl;

//            double chi = pow((experiment - theory),2)/pow(datapoint[i].error,2);
            return experiment - theory;

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

    // do the fit, obtain ra structure
    const APLCON::Result_t& ra = aplcon.DoFit();
//    cout << ra << endl;

    // Fits
    if (CEFTfit.GetNCalls() != 0)
    {
        cout << "Using CEFT:" << endl;
        cout << "Calls to Executable: " << CEFTfit.GetNCalls() << endl;
        cout << "Number of re-used folders: " << CEFTfit.GetNFolders() << endl;
    }
    else if (DISPfit.GetNCalls() != 0)
    {
        cout << "Using DISP:" << endl;
        cout << "Calls to Executable: " << DISPfit.GetNCalls() << endl;
        cout << "Number of re-used folders: " << DISPfit.GetNFolders() << endl;
    }

    // output time performance
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    cout << "*************************" << endl;
    cout << "Fit completed in " << elapsed_secs << " seconds." << endl;
    cout << "*************************" << endl;

    cout << "Instance " << instance << endl;
    cout << "E1E1 is found to be: " << ra.Variables.at("E1E1").Value.After << " +- " << ra.Variables.at("E1E1").Sigma.After << endl;

    cout << "M1M1 is found to be: " << ra.Variables.at("M1M1").Value.After << " +- " << ra.Variables.at("M1M1").Sigma.After << endl;

    cout << "E1M2 is found to be: " << ra.Variables.at("E1M2").Value.After << " +- " << ra.Variables.at("E1M2").Sigma.After << endl;

    cout << "M1E2 is found to be: " << ra.Variables.at("M1E2").Value.After << " +- " << ra.Variables.at("M1E2").Sigma.After << endl;

    cout << "Chi2 = " << ra.ChiSquare << endl;
    cout << "DOF = " << ra.NDoF << endl;
    cout << "Chi2/DOF = " << (ra.ChiSquare)/double(ra.NDoF) << endl;
    cout << "Probability = " << ra.Probability << endl;
    return ra;
}
