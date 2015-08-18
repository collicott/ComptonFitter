#include <iostream>
#include "CEFTwrapper.h"
#include "DISPwrapper.h"
#include "DataHandling.h"
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

struct set_params {
    int n_points;
    double e_min;
    double e_max;
    double t_min;
    double t_max;
    string data_type;

    set_params(int _n_points,
            double _e_min,
           double _e_max,
           double _t_min,
           double _t_max,
           string _data_type) :
        n_points(_n_points),
        e_min(_e_min),
        e_max(_e_max),
        t_min(_t_min),
        t_max(_t_max),
        data_type(_data_type)
    {}
};

class PseudoData {
public:
    PseudoData() {}

    vector<data> DataSet(set_params set);
    void Generate(vector<data> pseudo_data, params fitparam, bool smear, double dE, string theory_code);
};



int main(int argc, char *argv[])
{
    if (argc != 8)
    {
        cout << argc << endl;
        cout << "oh shit" << endl;
        return 1;
    }
    string theory_code = argv[1];

    cout << "Generating pseudo data: " << argv[2] << " "
                                       << argv[3] << " "
                                       << argv[4] << " "
                                       << argv[5] << " "
                                       << argv[6] << " "
                                       << argv[7] << endl;
    params fitparam(atof(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]));

    PseudoData pseudo;
//    params fitparam(11.2, 2.5, -3.3, 3.0, 0.2, 1.1);
//    params fitparam(11.2, 2.5, -1.0, 1.0, -2.5, 3.5);
    vector<set_params> set;

    set.push_back ({6,270.0,290.0,60,150,"Sigma_3"});
    set.push_back ({6,290.0,305.0,60,150,"Sigma_3"});
    set.push_back ({4,280.0,300.0,60,150,"Sigma_2x"});

    vector<data> pseudo_data;

//    for (auto i = 0; i < set.size(); i++)
//    {
//        vector<data> p = pseudo.DataSet(set[i]);
//        pseudo_data.insert(pseudo_data.end(), p.begin(), p.end());
//    }

    pseudo_data.push_back ({70, 275, 0, 0, "Sigma_3"});
    pseudo_data.push_back ({80, 275, 0, 0, "Sigma_3"});
    pseudo_data.push_back ({90, 275, 0, 0, "Sigma_3"});

    pseudo_data.push_back ({70, 295, 0, 0, "Sigma_3"});
    pseudo_data.push_back ({80, 295, 0, 0, "Sigma_3"});
    pseudo_data.push_back ({90, 295, 0, 0, "Sigma_3"});

    pseudo_data.push_back ({70, 295, 0, 0, "Sigma_2x"});
    pseudo_data.push_back ({80, 295, 0, 0, "Sigma_2x"});
    pseudo_data.push_back ({90, 295, 0, 0, "Sigma_2x"});


    //pseudo.Generate(pseudo_data, fitparam, true, 20.0, theory_code);
    if(theory_code == "all")
    {
        pseudo.Generate(pseudo_data, fitparam, false, 1.0, "Pascalutsa");
        pseudo.Generate(pseudo_data, fitparam, false, 1.0, "Pasquini");

    }
    else pseudo.Generate(pseudo_data, fitparam, false, 1.0, theory_code);

    return 0;

}

vector<data> PseudoData::DataSet(set_params set)
{
    // Set up random number distribution
    uniform_real_distribution<double> rdm_e(set.e_min, set.e_max); //(min, max)
    uniform_real_distribution<double> rdm_t(set.t_min, set.t_max);  //(min, max)

    mt19937 seed_e; seed_e.seed(random_device{}());
    mt19937 seed_t;  seed_t.seed(random_device{}());

    vector<data> DataSet;
    for (auto i = 0; i < set.n_points; i++)
    {
        double energy = rdm_e(seed_e);
        double theta  = rdm_t(seed_t);

        DataSet.push_back ({theta, energy, 0.0, 0.0, set.data_type});
    }

    return DataSet;
}

void PseudoData::Generate(vector<data> pseudo_data, params fitparam, bool smear, double percent_error, string theory_code)
{
    // Set up fitters
    CEFTwrapper CEFTfit;
    DISPwrapper DISPfit;



    // Generate pseudo
    for (auto i = 0; i < pseudo_data.size(); i++)
    {
         // Generate observable
         double theory = 0;
         if(theory_code == "Pascalutsa")
         theory = CEFTfit.Fit(pseudo_data[i].theta,
                              pseudo_data[i].energy,
                              fitparam.alpha,
                              fitparam.beta,
                              fitparam.E1E1,
                              fitparam.M1M1,
                              fitparam.E1M2,
                              fitparam.M1E2,
                              pseudo_data[i].data_type);
         else if (theory_code == "Pasquini")
         theory = DISPfit.Fit(pseudo_data[i].theta,
                              pseudo_data[i].energy,
                              fitparam.alpha,
                              fitparam.beta,
                              fitparam.E1E1,
                              fitparam.M1M1,
                              fitparam.E1M2,
                              fitparam.M1E2,
                              pseudo_data[i].data_type);

         if(smear)
         {
             // Set up random number distribution
             normal_distribution<double> rdm_t(theory, (percent_error/100.0)*theory);
             mt19937 seed_t; seed_t.seed(random_device{}());

             theory = rdm_t(seed_t);
         }


         // Generate error
         double error = theory*(percent_error/100.0);

         // Finally set pseudo data parameters
         pseudo_data[i].observable = theory;
         pseudo_data[i].error = error;
    }

    // Create file and store pseudo data
 /*   ofstream output;
    output.open("PseudoData.dat", std::ofstream::out);

    for (auto i = 0; i < pseudo_data.size(); i++)
    {

        output << pseudo_data[i].theta << " "
               << pseudo_data[i].energy << " "
               << pseudo_data[i].observable << " "
               << pseudo_data[i].error << " "
               << pseudo_data[i].data_type << " "
               << "lab" << endl;
    }*/

    FILE * output;
    string name = "data/PseudoData_" +  theory_code + ".dat";
    cout << "Pseudo data generated to file: " << name << endl;
    output = fopen (name.c_str(),"w");

    for (auto i = 0; i < pseudo_data.size(); i++)
    {
        fprintf(output, "%-8.2f %-8.2f %-8.6f %-8.6f %-0.10s lab\n",
                pseudo_data[i].theta,
                pseudo_data[i].energy,
                pseudo_data[i].observable,
                pseudo_data[i].error,
                pseudo_data[i].data_type.c_str());
    }

    fclose(output);

    return;
}
