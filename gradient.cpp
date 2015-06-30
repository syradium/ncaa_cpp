//
// Created by user on 6/30/15.
//

#include "gradient.h"

dvector GradE(dvector& thetaAngles, const dvector& phiAngles, dvector& magneticMoments, dvector& electronsNumber,
              const dvector& E0, const dvector& U0, const dmatrix& hopingIntegrals)
{
    cx_mat H(2 * L, 2 * L);
    H.fill(0);
    formHamiltonian(H, thetaAngles, phiAngles, electronsNumber, magneticMoments, E0, U0, hopingIntegrals, 0);

    //Make the matrix Hermitian
    //    H += H.t();
    //Finding the basis in which H is a diagonal transformation by calculating eigenvalues and eigenvectors
    cx_mat eigenVectors;
    vec eigenValues;
    eig_sym(eigenValues, eigenVectors, H);

    //2 * L means we have L atoms with two degrees of freedom each. First L terms terms theta and second L are phi polar angles
    dvector gradient(2 * L);

    for(unsigned int i = 0; i < gradient.size(); ++i)
    {
        //H derivative with respect to theta angles
        cx_mat dHdt(2 * L, 2 * L);
        //H derivative with respect to phi angles
        cx_mat dHdp(2 * L, 2 * L);

        dHdt.fill(0);
        dHdp.fill(0);

        double s1 = 0;
        double s2 = 0;

        //As mentioned [0; L) - are thetaAngles, [L; 2 * L) - phiAngles
        if(i < L)
        {
            dHdt(i, i)          = complex<double>( U0[i] * 0.5 * magneticMoments[i] * sin(thetaAngles[i]), 0);
            dHdt(i + L, i + L)  = complex<double>(-U0[i] * 0.5 * magneticMoments[i] * sin(thetaAngles[i]), 0);
            double x            = U0[i] * magneticMoments[i] * cos(thetaAngles[i]) * 0.5;
            dHdt(i, i + L)      = complex<double>(x * cos(phiAngles[i]), -x * sin(phiAngles[i]) );
            dHdt(i + L, i)      = complex<double>(x * cos(phiAngles[i]),  x * sin(phiAngles[i]) );

            //Moving dHdt to the basis of our Hamiltonian. H**T equals H inverse in case of such basis is a matrix of eigenvectors of a hermitian matrix is an unitary one
            cx_mat dHdt_newBasis = eigenVectors.t() * dHdt * eigenVectors;

            for(unsigned j = 0; j < 2 * L; ++j)
            {
                s1 += dHdt_newBasis(j, j).real();
                s2 += dHdt_newBasis(j, j).real() * atan( eigenValues(j) );
            }
        }
        else
        {
            //Introduce new shifted by L counter just for a convenience
            unsigned int j = i - L;

            double x            = U0[i] * magneticMoments[j] * sin(thetaAngles[j]) * 0.5;
            dHdp(j, j + L)      = complex<double>(-x * sin(phiAngles[j]), -x * cos(phiAngles[j]));
            dHdp(j + L, j)      = complex<double>(-x * sin(phiAngles[j]),  x * cos(phiAngles[j]));

            //Moving dHdp to the basis of our Hamiltonian. H**T equals H inverse in case of such basis
            cx_mat dHdp_newBasis = eigenVectors.t() * dHdp * eigenVectors;

            for(unsigned p = 0; p < 2 * L; ++p)
            {
                s1 += dHdp_newBasis(p, p).real();
                s2 += dHdp_newBasis(p, p).real() * atan( eigenValues(p) );
            }
        }
        //      Because of the lack of arcCtg function use trigonometrical equality between arctg and arcctg
        gradient[i] = (0.5 * s1 - s2 / boost::math::constants::pi<double>() );
    }
    return gradient;
}

