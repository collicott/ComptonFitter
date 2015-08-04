#include "DataHandling.h"
#include "CEFTwrapper.h"
#include "TLorentzVector.h"

#include <iostream>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <cstring>
using namespace std;

void DataHandling::process_file_list(const string &filename)
{
    ifstream input;

    input.open(filename);

    if(!input.is_open())
        throw result_file_error(filename);

    string newfile;
    while(input >> newfile)
    {
        process_file(newfile);
    }
    
    input.close();
}

void DataHandling::process_file(const string &filename)
{
//    cout << "processing file " << filename << endl;
    ifstream input;

    input.open(filename);

    if(!input.is_open())
        throw result_file_error(filename);

    double theta;
    double energy;
    double observable;
    double error;
    string observable_type;
    string frame;
    while(input >> theta >> energy >> observable >> error >> observable_type >> frame)
    {
 //       cout << theta <<" "<< energy <<" "<<  observable <<" "<<  error <<" "<<  observable_type << " " << frame;

        // Use toupper to avoid case sensitivity
        std::transform(observable_type.begin(),observable_type.end(),observable_type.begin(), ::toupper);
        std::transform(frame.begin(),frame.end(),frame.begin(), ::toupper);

        // if CM frame data point, convert to LAB
        if (frame == "CM")
        {
            theta = Convert_CM_to_LAB(theta,energy);
        }
//        else cout << endl;

        if (observable_type == "SIGMA_3")
            datalist.push_back ({theta, energy, observable, error, "Sigma_3"});

        else if (observable_type == "SIGMA_2X")
            datalist.push_back ({theta, energy, observable, error, "Sigma_2x"});

        else if (observable_type == "SIGMA_2z")
            datalist.push_back ({theta, energy, observable, error, "Sigma_2z"});

        else if (observable_type == "CROSS")
            datalist.push_back ({theta, energy, observable, error, "Cross"});

        else cout << "Warning: skipping *unknown* data type " << observable_type << endl;
    }
    input.close();
}

double DataHandling::Convert_CM_to_LAB(double th_lab, double E_beam)
{
    beam   = TLorentzVector(0.0, 0.0 , E_beam, E_beam );
    target = TLorentzVector(0.0, 0.0 , 0.0, 938.27 );

    particle = TLorentzVector(0.0, 0.0, 1.0, 1.0);
    particle.SetTheta(th_lab*TMath::DegToRad());

    TVector3 BoostVect;
    //use -1.0 for LAB to CM
    BoostVect = 1.0 * (beam + target).BoostVector();
    particle.Boost(BoostVect);

    return particle.Theta()*TMath::RadToDeg();
}

