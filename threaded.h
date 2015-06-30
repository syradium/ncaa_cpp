//
// Created by user on 6/30/15.
//

#ifndef DIPLOMA_THREADED_H
#define DIPLOMA_THREADED_H

#include "definitions.h"

void SConsistThreaded(int i, std::pair< dvector, dvector >& result, const std::pair< dvector, dvector >& angles, const dvector& M,
                      const dvector& N, const std::pair< dvector, dvector >& E0U0, const dmatrix& hopingIntegrals, dvector& E, bool& isConsist);
#endif //DIPLOMA_THREADED_H
