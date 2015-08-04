#include "DisplayResults.h"

#include <iostream>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <cstring>
using namespace std;

void DisplayResults::ShowFitResults(vector<data> datapoint, const double a, const double b, const double E1E1, const double M1M1, const double E1M2, const double M1E2)
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

    for (auto i = 0; i < energy_list.size(); i++)
    {
        cout << energy_list[i] << " " << observable_list[i] << endl;
    }

}
