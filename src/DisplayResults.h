#ifndef DisplayResults_H
#define DisplayResults_H

#include "CEFTwrapper.h"
#include "DISPwrapper.h"
#include "DataHandling.h"
#include "TFile.h"

#include <string>
#include <fstream>
#include <stdlib.h>
#include <exception>
#include <vector>

class DisplayResults {
private:

public:

    DisplayResults() {}

    void ShowFitResults(std::vector<data> datapoint, const double a, const double b, const double E1E1, const double M1M1, const double E1M2, const double M1E2);
};

#endif
