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

// Test phils fit - datapoints not fixed
int main(int argc, char *argv[])
{
    // Set up some time
    clock_t begin = clock();

    // Quickly parse through arguments to find theory and filelist
    string theory_code;
    string filelist;

    string flag(argv[2]);
    if (flag == "--Pascalutsa") theory_code = "Pascalutsa";
    else if (flag == "--Pasquini") theory_code = "Pasquini";
    else {
        cout << "Error: unknown theory " << argv[2] << endl;
        return 1;
    }

    // Set up Data handler
    DataHandling dataHandler;
    dataHandler.process_file_list(argv[1]);
    vector<data> datapoint = dataHandler.GetDataList();

    // Set up fitters
    CEFTwrapper CEFTfit;
    DISPwrapper DISPfit;

    // Set up some ROOT storage
    TFile *f = new TFile("phil_fit.root","RECREATE");
    TNtuple *ntuple = new TNtuple("phil_fit",
                                  "phil_fit",
                                  "call:alpha:beta:E1E1:M1M1:E1M2:M1E2");

    int call = 0;

    // Set up APLCON
    APLCON::Fit_Settings_t settings = APLCON::Fit_Settings_t::Default;
    settings.MaxIterations = 1000;
    //settings.UnmeasuredStepSizeFactor = 1;
    settings.DebugLevel = 5;
    APLCON aplcon("Polarisabilities", settings);

    // Set initial values

    // old alpha-beta
//    params fitparam(12.1, 1.6, -4.3, 2.9, -0.02, 2.2);   // (nominal HDPV) -3.7
//    params fitparam(12.1, 1.6, -3.3, 3.0, 0.2, 1.1);     // (nominal BxPT) -3.7

    // new alpha-beta
//      params fitparam(11.2, 2.5, -4.3, 2.9, -0.02, 2.2);    // (nominal HDPV)   -3.6
    params fitparam(11.2, 2.5, -3.3, 3.0, 0.2, 1.1);     // (nominal BxPT) -3.7

    // Add params as unmeasured variables
    aplcon.AddUnmeasuredVariable("E1E1",fitparam.E1E1);

    // Add data points as measured variables
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
        auto equality_constraint = [&CEFTfit,&DISPfit,&fitparam,&datapoint,i,&ntuple,&call,&theory_code]
                (double experiment,
                 double E1E1)
        {
            double theory = 0;
            if(theory_code == "Pascalutsa")
                 theory = CEFTfit.Fit(datapoint[i].theta,
                                      datapoint[i].energy,
                                      fitparam.alpha,
                                      fitparam.beta,
                                      E1E1,
                                      fitparam.M1M1,
                                      fitparam.E1M2,
                                      fitparam.M1E2,
                                      datapoint[i].data_type);
            else if (theory_code == "Pasquini")
                 theory = DISPfit.Fit(datapoint[i].theta,
                                      datapoint[i].energy,
                                      fitparam.alpha,
                                      fitparam.beta,
                                      E1E1,
                                      fitparam.M1M1,
                                      fitparam.E1M2,
                                      fitparam.M1E2,
                                      datapoint[i].data_type);
            ntuple->Fill(call,
                         fitparam.alpha,
                         fitparam.beta,
                         E1E1,
                         fitparam.M1M1,
                         fitparam.E1M2,
                         fitparam.M1E2);
            call++;

//          double chi = pow((experiment - theory),2)/pow(datapoint[i].error,2);
            return experiment - theory;
        };

        aplcon.AddConstraint(
                    constraintname.str(),
                    {pointname.str(),"E1E1"},
                     equality_constraint);
    }

    // do the fit, obtain ra structure
    const APLCON::Result_t& ra = aplcon.DoFit();
    cout << ra << endl;

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
    cout << "E1E1 is found to be: " << ra.Variables.at("E1E1").Value.After << " +- " << ra.Variables.at("E1E1").Sigma.After << endl;


    // Let's display the results
    DisplayResults display;
    display.ShowFitResults(f,
                           datapoint,
                           fitparam.alpha,
                           fitparam.beta,
                           ra.Variables.at("E1E1").Value.After,
                           fitparam.M1M1,
                           fitparam.E1M2,
                           fitparam.M1E2,
                           theory_code);
    f->Write();

    return 0;
}
