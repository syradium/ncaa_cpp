//
// Created by user on 6/30/15.
//

#ifndef DIPLOMA_DEFINITIONS_H
#define DIPLOMA_DEFINITIONS_H
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/thread/thread.hpp>
#include <boost/numeric/ublas/hermitian.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/foreach.hpp>
#include <boost/math/constants/constants.hpp>
#include <iostream>
#include <complex>
#include <vector>
#include <iomanip>
#include <limits>
#include <fstream>
#include <sstream>
#include <numeric>
#include <armadillo>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <vector>

#include "constants.h"

using std::complex;
using std::cout;
using std::numeric_limits;
using namespace boost::numeric::ublas;
using namespace boost::math::tools;
using namespace arma;

typedef std::vector<std::vector<double> > dmatrix;
typedef std::vector< double > dvector;
typedef hermitian_matrix<std::complex<double>, upper> herm_matrix;

//H - hamiltonian, tAngle - thetta angles, pAngle - phi angles, N - number of d-electrons on each atom, M - magnetic moments of each atom
template <typename T> void formHamiltonian(T& H, const dvector& tAngle, const dvector& pAngle, const dvector& N, const dvector& M, const dvector& E0, const dvector& U0, const dmatrix& V, const double& dh) {
    //hermitian_matrix<std::complex<double>, upper> H (2 * L, 2 * L);
    for (int i = 0; i < L; ++i) {
        H(i, i)         = complex<double>(E0[i] - dh + U0[i] * 0.5 * (N[i] - M[i] * cos(tAngle[i])), 0);
        H(i + L, i + L) = complex<double>(E0[i] + dh + U0[i] * 0.5 * (N[i] + M[i] * cos(tAngle[i])), 0);
        for (int j = i + 1; j < L; ++j)
            H(i + L, j + L) = H(i, j) = complex<double>(V[i][j], 0);
        double x = U0[i] * M[i] * sin(tAngle[i]) * 0.5;
        H(i, i + L) = complex<double>(x * cos(pAngle[i]), -x * sin(pAngle[i]));
    }
}

matrix<double> Green(const dvector& a, const dvector& b);
double CN(const matrix<double>& t);
double CalculateEnergy(const matrix<double>& t);
std::pair< dvector, dvector > ThreeDiag(const hermitian_matrix<std::complex<double>, upper>& H, int pIndex, int qIndex, const complex<double>& k, const complex<double>& m);

#endif //DIPLOMA_DEFINITIONS_H
