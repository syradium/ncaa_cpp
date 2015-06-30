//
// Created by user on 6/30/15.
//

#ifndef DIPLOMA_GRADIENT_H
#define DIPLOMA_GRADIENT_H

#include "definitions.h"
#include "constants.h"

typedef std::vector< double > dvector;
typedef std::vector<std::vector<double> > dmatrix;
using namespace arma;
using std::complex;


dvector GradE(dvector& thetaAngles, const dvector& phiAngles, dvector& magneticMoments, dvector& electronsNumber,
              const dvector& E0, const dvector& U0, const dmatrix& hopingIntegrals);


#endif //DIPLOMA_GRADIENT_H
