#ifndef CEFTWRAPPER_H
#define CEFTWRAPPER_H

#include <string>
#include <fstream>
#include <stdlib.h>
#include <exception>
#include <sys/stat.h>
#include <map>

class Fitter {
private:

    unsigned int n_calls = 0;
    unsigned int n_folders = 0;

    // results map
    std::map<std::string, double> results_map;

    // The command to call
    std::string folder;
    std::string command = "sh ComptonEFT/runComptonEFT";


public:

    Fitter() {}

    double Fit(const double th, const double E, const double a, const double b, const double E1E1, const double M1M1, const double E1M2, const double M1E2, std::string data_type);

    unsigned int GetNCalls() const { return n_calls; }
    unsigned int GetNFolders() const { return n_folders; }

    double read_double(const std::string& filename) const;
    double GetSigma2x() const;
    double GetSigma2z() const;
    double GetSigma3() const;
    double GetCross() const;

    // Error handling
    class result_file_error: public std::exception {
    protected:
        std::string filename;
    public:
        result_file_error(const std::string& file): filename(file) {}
    public:
        const char *what() const throw() { return filename.c_str(); }
    };

    class fit_command_error: public std::exception {
    protected:
        std::string msg;
        int code;
    public:
        fit_command_error(const int c): code(c), msg("Fit command return code "+std::to_string(c)) {}
    public:
        const char *what() const throw() { return msg.c_str(); }
        const int ErrorCode() const throw() { return code; }
    };
};

#endif
