#include <iostream>
#include <APLCON.hpp>
#include "CEFTwrapper.h"
#include "DataHandling.h"
#include "TFile.h"
#include "TNtuple.h"
#include <ctime>

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

     // Set up fitter
    Fitter myfit;

    // Set initial values

    // old alpha-beta
//      params fitparam(12.1, 1.6, -3.3, 3.0, 0.2, 1.1);     // (nominal BxPT) -3.7
//    params fitparam(12.1, 1.6, -4.3, 2.9, -0.02, 2.2);   // (nominal HDPV) -3.7

    // new alpha-beta
//    params fitparam(11.2, 2.5, -3.3, 3.0, 0.2, 1.1);     // (nominal BxPT) -3.7
    params fitparam(11.2, 2.5, -4.3, 2.9, -0.02, 2.2);    // (nominal HDPV)   -3.6

    // Add params as unmeasured variables
    aplcon.AddUnmeasuredVariable("E1E1",fitparam.E1E1);

    // Set up Data handler
    DataHandling dataHandler;

    string filelist(argv[1]);
    dataHandler.process_file_list(filelist);
    vector<data> datapoint = dataHandler.GetDataList();

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
        auto equality_constraint = [&myfit,&fitparam,&datapoint,i,&ntuple,&call]
                (double experiment,
                 double E1E1)
        {
            auto theory = myfit.Fit(datapoint[i].theta,
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

    cout << "Calls to Executable: " << myfit.GetNCalls() << endl;
    cout << "Number of re-used folders: " << myfit.GetNFolders() << endl;

    f->Write();

    // output time performance
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    cout << "*************************" << endl;
    cout << "Fit completed in " << elapsed_secs << " seconds." << endl;

    return 0;
}
