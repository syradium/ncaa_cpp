//
// Created by user on 6/29/15.
//

#ifndef DIPLOMA_CONSTANTS_H
#define DIPLOMA_CONSTANTS_H

const int CORE_COUNT = 1;
const int L = 3;
//The maximum number of iteration in selfconsist procedure. If that number is reached we consider that selfconsist procedure will never end, so we terminate it.
const int maxIterCount = 400;
const double delta  = 1e-10;
const double eps    = 1e-11 + std::numeric_limits<double>::epsilon();

#endif //DIPLOMA_CONSTANTS_H
