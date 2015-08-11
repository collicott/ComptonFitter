#include "DisplayResults.h"

#include <iostream>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <cstring>
using namespace std;

void DisplayResults::ShowFitResults(TFile* f, vector<data> datapoint, const double alpha, const double beta, const double E1E1, const double M1M1, const double E1M2, const double M1E2, string theory_code)
{
    // Let's grab a new instance of each fitter
    CEFTwrapper CEFTfit;
    DISPwrapper DISPfit;

    // Let's run through the datalist once to get a set of energies/observables
    vector<double> energy_list;
    vector<string> observable_list;
    for (auto i = 0; i < datapoint.size(); i++)
    {
        // New energy? Store energy and data type
        if(find(energy_list.begin(), energy_list.end(), datapoint[i].energy) == energy_list.end())
        {
            energy_list.push_back(datapoint[i].energy);
            observable_list.push_back(datapoint[i].data_type);
        }
    }

    double chi_sum = 0.0;
    for (auto i = 0; i < energy_list.size(); i++)
    {
        cout << energy_list[i] << " " << observable_list[i] << endl;
        vector<double> experimental_points;
        vector<double> experimental_errors;
        vector<double> theory_points;
        vector<double> theta_points;
        vector<double> theta_errors;

        for (auto j = 0; j < datapoint.size(); j++)
        {
            if (datapoint[j].energy != energy_list[i]) continue;

            double theory = 0;
            if(theory_code == "Pascalutsa")
                theory = CEFTfit.Fit(datapoint[j].theta,
                                     datapoint[j].energy,
                                     alpha,
                                     beta,
                                     E1E1,
                                     M1M1,
                                     E1M2,
                                     M1E2,
                                     datapoint[j].data_type);
            else if (theory_code == "Pasquini")
                theory = DISPfit.Fit(datapoint[j].theta,
                                     datapoint[j].energy,
                                     alpha,
                                     beta,
                                     E1E1,
                                     M1M1,
                                     E1M2,
                                     M1E2,
                                     datapoint[j].data_type);

            double chi = pow((datapoint[j].observable - theory),2)/pow(datapoint[j].error,2);
            chi_sum+=chi;

            cout << energy_list[i] << ": " << datapoint[j].data_type << " " << datapoint[j].observable << " -- Theory: " << theory << " -- Chi: " << chi << endl;


            // Set up the TGraphs
            theta_points.push_back(datapoint[j].theta);
            theta_errors.push_back(0.0);
            experimental_points.push_back(datapoint[j].observable);
            experimental_errors.push_back(datapoint[j].error);
            theory_points.push_back(theory);

        }

        // Make dem TGraphs
        string name = observable_list[i] + "_" + to_string(energy_list[i]);
        TCanvas *c1 = new TCanvas(name.c_str(),name.c_str(),200,10,700,500);

        TGraphErrors *exp = new TGraphErrors(theta_points.size(), &theta_points[0], &experimental_points[0],&theta_errors[0],&experimental_errors[0]);
        exp->SetMarkerStyle(21);
        exp->SetMarkerColor(kBlue);
        exp->SetLineColor(kBlue);
        exp->SetTitle(name.c_str());
        exp->GetXaxis()->SetLimits(0.0,180.0);

        TGraph *thr = new TGraph(theta_points.size(), &theta_points[0], &theory_points[0]);
        thr->SetMarkerStyle(3);
        thr->SetMarkerColor(kRed);

        exp->Draw("AP");
        thr->Draw("P");

        c1->Write();

    }
    struct constraints {
        double value;
        double error;

        constraints(double _value,
               double _error) :
            value(_value),
            error(_error)
        {}
    };

    // Add in some chi punishment for the four constraints
    double chi;
//    constraints alpha_beta_diff(10.5, 1.6);  // Old numbers
    constraints alpha_beta_diff(7.6, 1.7); // Griesshammer
    constraints alpha_beta_sum(13.8, 0.4);
    constraints gamma_0(-1.01, 0.08);
    constraints gamma_pi(8.0, 1.8);

    cout << endl;
    chi = pow((alpha + beta - alpha_beta_sum.value),2)/pow(alpha_beta_sum.error,2);
    chi_sum+=chi;
    cout << "ab_sum: " << alpha + beta << " --Chi: " << chi << endl;

    chi = pow((alpha - beta - alpha_beta_diff.value),2)/pow(alpha_beta_diff.error,2);
    chi_sum+=chi;
    cout << "ab_diff: " << alpha - beta << " --Chi: " << chi << endl;

    chi = pow((-E1E1 - M1M1 - E1M2 - M1E2 - gamma_0.value),2)/pow(gamma_0.error,2);
    chi_sum+=chi;
    cout << "gamma_0: " << -E1E1 - M1M1 - E1M2 - M1E2 << " --Chi: " << chi << endl;

    chi = pow((-E1E1 + M1M1 - E1M2 + M1E2 - gamma_pi.value),2)/pow(gamma_pi.error,2);
    chi_sum+=chi;
    cout << "gamma_pi: " << -E1E1 + M1M1 - E1M2 + M1E2  << " --Chi: " << chi << endl;

    cout << "******************************************" << endl;
    cout << "Chi_squared = "<< chi_sum << endl;

}
