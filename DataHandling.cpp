#include "DataHandling.h"
#include "CEFTwrapper.h"

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
    cout << "processing file " << filename << endl;
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
        cout << theta <<" "<< energy <<" "<<  observable <<" "<<  error <<" "<<  observable_type << " " << frame << endl;

        // Use toupper to avoid case sensitivity
        std::transform(observable_type.begin(),observable_type.end(),observable_type.begin(), ::toupper);
        std::transform(frame.begin(),frame.end(),frame.begin(), ::toupper);

        // if LAB frame, convert to CM for Pascalutsa
   //     if (frame == "LAB")
     //       theta = Convert_LAB_to_CM(theta);

        if (observable_type == "SIGMA_3")
            datalist.push_back ({theta, energy, observable, error, Sigma_3});

        else if (observable_type == "SIGMA_2X")
            datalist.push_back ({theta, energy, observable, error, Sigma_2x});

        else if (observable_type == "SIGMA_2z")
            datalist.push_back ({theta, energy, observable, error, Sigma_2z});

        else cout << "Warning: skipping *unknown* data type " << observable_type << endl;
    }
    input.close();
}

/*
vector<data> DataHandling::GetDataList()
{
    vector<data> datalist;

    datalist.push_back ({75, 297, 0.2697, 0.0388, Sigma_3});
    datalist.push_back ({85, 297, 0.2551, 0.0274, Sigma_3});
    datalist.push_back ({95, 297, 0.2228, 0.267, Sigma_3});
    datalist.push_back ({107.5, 297, 0.1597, 0.0316, Sigma_3});

    cout << datalist[0].energy << endl;

    return datalist;
}
*/

/*double DataHandling::Convert_LAB_to_CM(double th_lab)
{
    return th_lab;
}
*/
