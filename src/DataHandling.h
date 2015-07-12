#ifndef DATAHANDLING_H
#define DATAHANDLING_H

#include "CEFTwrapper.h"
#include "TLorentzVector.h"

#include <string>
#include <fstream>
#include <stdlib.h>
#include <exception>
#include <vector>

struct data {
    double theta;
    double energy;
    double observable;
    double error;
    int    data_type;


    data(double _theta,
         double _energy,
         double _observable,
         double _error,
         int    _data_type) :
        theta(_theta),
        energy(_energy),
        observable(_observable),
        error(_error),
        data_type(_data_type)
    {}
};

enum
{
    Sigma_2x = 1,
    Sigma_2z = 2,
    Sigma_3  = 3,
    Cross    = 4,
};

class DataHandling {
private:

    std::vector<data> datalist;

    // CM conversion stuff
    TLorentzVector beam;
    TLorentzVector target;
    TLorentzVector particle;


public:

    class result_file_error: public std::exception {
    protected:
        std::string filename;
    public:
        result_file_error(const std::string& file): filename(file) {}

    public:
        const char *Data_Handling() const throw() { return filename.c_str(); }
    };

    DataHandling() {}

    const std::vector<data> GetDataList() {return datalist;}
    void process_file_list(const std::string& filename);
    void process_file(const std::string& filename);
    double Convert_CM_to_LAB(double th_lab, double E_beam);
};

#endif
