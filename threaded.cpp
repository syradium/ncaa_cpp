//
// Created by user on 6/30/15.
//

#include "threaded.h"
bool interrupt = false;

void SConsistThreaded(int i, std::pair< dvector, dvector >& result, const std::pair< dvector, dvector >& angles, const dvector& M,
                      const dvector& N, const std::pair< dvector, dvector >& E0U0, const dmatrix& hopingIntegrals, dvector& E, bool& isConsist)
{
    double magneticField = 0.0;
    //Saving magnetic moments for each atom
    dvector M1(M); dvector M0(M);
    //Saving number of d-electrons on each atom
    dvector N1(N); dvector N0(N);
    //Hamiltonian
    hermitian_matrix< complex<double>, upper> H(2 * L, 2 * L);
    //The matrix of P and Q coefficients of the fraction that represents green`s function`s matrix element on i-th atom
    matrix<double> greenFraction;
    //    Number of iterations used
    dvector E0 = E0U0.first;
    dvector U0 = E0U0.second;
    unsigned int itr = 0;

    do
    {
        ++itr;
        //Copy new values of the number d-electrons and magnetic moments and proceed the process
        N0[i] = N1[i];
        M0[i] = M1[i];

        H.clear(); greenFraction.clear();

        formHamiltonian(H, angles.first, angles.second, N0, M0, E0, U0, hopingIntegrals, magneticField);

        std::pair< dvector, dvector > cofficients = ThreeDiag(H, i, i, complex<double>(1, 0), complex<double>(1, 0));

        int ii0 = 0;
        greenFraction = Green(cofficients.first, cofficients.second, ii0);
        double Nu = CN(greenFraction, cofficients.first.size() - 1);
        double Eu = CalculateEnergy(greenFraction, cofficients.first.size() - 1);

        cofficients = ThreeDiag(H, i + L, i + L, complex<double>(1, 0), complex<double>(1, 0));
        greenFraction = Green(cofficients.first, cofficients.second, ii0);

        double Nd = CN(greenFraction, cofficients.first.size() - 1);
        double Ed = CalculateEnergy(greenFraction, cofficients.first.size() - 1);

        cofficients = ThreeDiag(H, i, i + L, complex<double>(1, 0), complex<double>(1, 0));

        greenFraction = Green(cofficients.first, cofficients.second, ii0);

        double SFp = CN(greenFraction, cofficients.first.size() - 1);

        cofficients = ThreeDiag(H, i, i + L, complex<double>(1, 0), complex<double>(-1, 0));
        greenFraction = Green(cofficients.first, cofficients.second, ii0);

        double SFn = CN(greenFraction, cofficients.first.size() - 1);

        cofficients = ThreeDiag(H, i, i + L, complex<double>(0, 1), complex<double>(1, 0));
        greenFraction = Green(cofficients.first, cofficients.second, ii0);
        double AFp = CN(greenFraction, cofficients.first.size() - 1);

        cofficients = ThreeDiag(H, i, i + L, complex<double>(0, 1), complex<double>(-1, 0));
        greenFraction = Green(cofficients.first, cofficients.second, ii0);
        double AFn = CN(greenFraction, cofficients.first.size() - 1);

        N1[i] = Nu + Nd;
        M1[i] = (Nu - Nd) * cos(angles.first[i]) - ((SFp - SFn) * cos(angles.second[i]) - (AFp - AFn) * sin(angles.second[i])) * sin(angles.first[i]);
        E[i] = Eu + Ed - U0[i] * ( N1[i] * N1[i] - M1[i] * M1[i] ) * 0.25 ;
    }
    while ( (std::abs(M1[i] - M0[i]) > delta || std::abs(N1[i] - N0[i]) > delta) && itr < maxIterCount ) ;
    result.first[i] = N1[i];
    result.second[i] = M1[i];

    if( itr == maxIterCount )
    {
        interrupt = true;
        return;
    }
    isConsist &= (itr == 1);
}