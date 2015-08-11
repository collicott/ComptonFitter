#ifndef DisplayResults_H
#define DisplayResults_H

#include "CEFTwrapper.h"
#include "DISPwrapper.h"
#include "DataHandling.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TAxis.h"

#include <string>
#include <fstream>
#include <stdlib.h>
#include <exception>
#include <vector>

class DisplayResults {
private:

public:

    DisplayResults() {}

    void ShowFitResults(TFile* f, std::vector<data> datapoint, const double alpha, const double beta, const double E1E1, const double M1M1, const double E1M2, const double M1E2, std::string theory_code);
};

#endif
