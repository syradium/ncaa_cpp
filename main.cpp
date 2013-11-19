/* 
 * File:   main.cpp
 * Author: Radium
 *
 * Created on January 14, 2013, 1:15 AM
 */

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
#include <c++/4.7/complex>

using namespace arma;


const int CORE_COUNT = 4;

bool interrupt = false;

const int L = 3;
//The maximum number of iteration in selfconsist procedure. If that number is reached we consider that selfconsist procedure will never end, so we terminate it.
const int maxIterCount = 400;

using std::complex;
using std::numeric_limits;
using namespace boost::numeric::ublas;
using namespace boost::math::tools;

typedef std::vector<std::vector<double> > dmatrix;
typedef std::vector< double > dvector;
typedef hermitian_matrix<std::complex<double>, upper> herm_matrix;

const double delta  = 1e-10;
const double eps    = 1e-11 + std::numeric_limits<double>::epsilon(); 

#include "debug_functions.h"


dvector GradE(dvector& thetaAngles, const dvector& phiAngles, dvector& magneticMoments, dvector& electronsNumber, const dvector& E0, const dvector& U0, const dmatrix& hopingIntegrals);

//H - hamiltonian, tAngle - thetta angles, pAngle - phi angles, N - number of d-electrons on each atom, M - magnetic moments of each atom
template <typename T> void formHamiltonian(T& H, const dvector& tAngle, const dvector& pAngle, const dvector& N, const dvector& M, const dvector& E0, const dvector& U0, const dmatrix& V, const double& dh)
{
    //hermitian_matrix<std::complex<double>, upper> H (2 * L, 2 * L);
    for (int i = 0; i < L; ++i)
    {
        H(i, i)         = complex<double>(E0[i] - dh + U0[i] * 0.5 * (N[i] - M[i] * cos(tAngle[i])), 0);
        H(i + L, i + L) = complex<double>(E0[i] + dh + U0[i] * 0.5 * (N[i] + M[i] * cos(tAngle[i])), 0);
        for (int j = i + 1; j < L; ++j)
            H(i + L, j + L) = H(i, j) = complex<double>(V[i][j], 0);
        double x = U0[i] * M[i] * sin(tAngle[i]) * 0.5;
        H(i, i + L) = complex<double>(x * cos(pAngle[i]), -x * sin(pAngle[i]));        
    }
}

std::pair< dvector, dvector > ThreeDiag(const hermitian_matrix<std::complex<double>, upper>& H, int pIndex, int qIndex, const complex<double>& k, const complex<double>& m)
{
    vector< std::complex<double> > y0(2 * L, 0);
    vector< std::complex<double> > y1(2 * L, 0);
    vector< std::complex<double> > y(2 * L, 0);
    //Subdiagonal coefficients of matrix
    std::vector< std::complex<double> > a(2 * L, 0);
    //Diagonal coefficiencrs of matrix
    std::vector< std::complex<double> > b(2 * L - 1, 0);
    /*
                                                  p              q
            First basis vector. y0 = ( 0, ..., 0, k, 0, ...., 0, m, 0 ... ) 
     */
    y0[pIndex] = k;
    y0[qIndex] = m;
    y0 /= norm_2(y0);
    //Acting on y0 by H - ( H * y0 )
    vector< std::complex<double> > Hy = prod(H, y0);
    //First coefficient a of threediagonalized matrix ( < y0,  (H * y0) > ) 
    a[0] = inner_prod(conj(y0), Hy);
    vector< std::complex<double> > buf = a[0] * y0;
    //Caution: y is not actually a basis vector. Here it is a y1 * b
    y = Hy - buf;
    double x = norm_2(y);
    int imax;
    for (imax = 0; imax < (2 * L - 1) && x * x > /*1e-8*/ ::eps * 1e4; ++imax)
    {
        b[imax] = x;
        //Now it is a basis
        y1 = y / b[imax];
        Hy = prod(H, y1);
        a[imax + 1] = inner_prod(conj(y1), Hy);
        y = Hy - b[imax] * y0 - a[imax + 1] * y1;
        y0 = y1;
        x = norm_2(y);
    }
    std::vector< double > aa;
    std::vector< double > bb;

    BOOST_FOREACH(std::complex<double> i, a)
    {
        aa.push_back(i.real());
    }
    BOOST_FOREACH(std::complex<double> i, b)
    {
        bb.push_back(i.real() * i.real() );
    }
    aa.resize(imax + 1);
    return std::make_pair(aa, bb);
}

double f(const double& x, const double& a0, const double& b0, const std::vector<double>& p0, const std::vector<double>& q0, const int& ii0)
{
    double s = 0;
    for (int i = 0; i <= ii0/*&& std::abs(p0[i]) >= 1E-15*/ /*0.00000001*/; ++i)
//        if( std::abs(p0[i]) >= 1E-15 )
                s += p0[i] / (x - q0[i]);
    return x - a0 - b0 * s;
}

double bisection(const double& a, const double& b, const double& a0, const double& b0, const std::vector<double>& p0, const std::vector<double>& q0, const int& ii0)
{

    double l = a;
    double r = b;
    unsigned iter = 0;
    while (r - l > ::eps )
    {
        double m = (l + r) / 2;
        if(++iter > 50)
                return a;

        if( f(m, a0, b0, p0, q0, ii0) * f(r, a0, b0, p0, q0, ii0) <= 0 )
            l = m;
        else
            r = m;
    }
    return l;

}

void FR(const double& a, const double& b, matrix<double>& t, std::vector<double>& p0, std::vector<double>& q0, int& ii0)
{

    double a0 = a;
    double b0 = b; //* b;
    std::vector<double> z(t.size2(), 0);

    int firstInfinitesimal = 0;
    while( firstInfinitesimal < 2 * L && std::abs(t(0, firstInfinitesimal) * b0) <= ::eps * 1e4 )
        ++firstInfinitesimal;

    if( firstInfinitesimal == 2 * L )
    {
            t(0, 0) = 1;
            t(1, 0) = a;
            ii0 = 0;
            return;
    }
    
    p0[0] = t(0, firstInfinitesimal);
    q0[0] = t(1, firstInfinitesimal);

    int j = 0;
    for(int i = firstInfinitesimal + 1; i <= ii0; ++i )
    {
        if( std::abs(t(0, i) * b0) >= ::eps * 1e4)
        {
            if(std::abs(t(1, i) - q0[j]) >= ::eps * 1e4)
            {
                ++j;
                p0[j] = t(0, i);
                q0[j] = t(1, i);
            }
            else   
                p0[j] = p0[j] + t(0, i);
        }
    }

    ii0 = j;
    
    double leftBound = q0[0]- 1E4;
    double rightBound = q0[ii0] + 1E4;
    
    z[0] = bisection( leftBound, q0[0] - ::eps * 10, a0, b0, p0, q0, ii0);

    for (int i = 0; i <= ii0; ++i)
        z[i + 1] = bisection(q0[i] + ::eps * 10, q0[i + 1] - ::eps * 10, a0, b0, p0, q0, ii0);
    z[ii0 + 1] = bisection(q0[ii0] + ::eps * 10, rightBound, a0, b0, p0, q0, ii0);

    for (int i = 0; i <= ii0 + 1; ++i)
    {
        double pTerm = 1;
        double pDenom = 1;
        for (int j = 0; j <= ii0; ++j)
        {
            pTerm *= z[i] - q0[j];
            if (i != j)
                pDenom *= z[i] - z[j];
        }
        if (i != ii0 + 1)
            pDenom *= z[i] - z[ii0 + 1];
        t(0, i) = pTerm / pDenom;
        t(1, i) = z[i];
    }    
     ++ii0;   
}


//Calculates the coefficients P and Q of the Green`s function`s matrix element represented by a fraction
//"a" and "b" are the diagonal and subdiagonal elemnts of the hamiltonian

matrix<double> Green(const std::vector<double>& a, const std::vector<double>& b, int& ii0)
{
    int imax = a.size() - 1;
    matrix<double> t(2, 2 * L);
    std::vector<double> p0(2 * L, 0);
    std::vector<double> q0(2 * L, 0);

    for (int i = 0; i < 2 * L; ++i)
        t(1, i) = t(0, i) = 0;
    t(0, 0) = 1;
    //t(1,0) = a[imax];		
    t(1, 0) = a.back();

    int ii00 = 0;
    for (int i = imax - 1; i >= 0; --i)
        FR(a[i], b[i], t, p0, q0, ii00);
    ii0 = ii00;
    return t;
}

//Carry out an integration of the density matrix, which is the number of d-electrons, if a diagonal element is considered
double CN(const matrix<double>& t, const int& imax)
{
    double s1 = 0;
    double s2 = 0;
    for (int i = 0; i <= imax; ++i)
    {
        s1 += t(0, i);
        s2 += t(0, i) * std::atan(t(1, i));
    }
    
//    std::cout << "sum(s1): " << s1 << std::endl;
    return 0.5 - s2 * 0.318309886183790671;
}

double CalculateEnergy(const matrix<double>& t, const int& imax)
{
    double s1 = 0;
    double s2 = 0;
    double s3 = 0;
    
    for (int i = 0; i <= imax; ++i)
    {
        double PQ = t(0, i) * t(1, i);
        s1 += PQ;
        s2 += PQ * atan( t(1, i) );
        s3 += t(0, i) * log(t(1, i) * t(1, i) + 1);
    }
    return 0.5 * s1 - (s2 - 0.5 * s3) * 0.318309886183790671;
}

void GetParam(std::ifstream& file, std::vector<double>& arr)
{
    std::string buf;
    std::getline(file, buf);
    std::istringstream line(buf);
    std::copy(std::istream_iterator<double>(line), std::istream_iterator<double>(), std::back_inserter(arr));
}

bool SConsist(int i, const dvector& tAngle, const dvector& pAngle, dvector& M, dvector& N, const dvector& E0, const dvector& U0, const dmatrix& hopingIntegrals, dvector& E)
{
    double magneticField = 0;
    //Saving magnetic moments for each atom
    dvector M1(M);
    //Saving number of d-electrons on each atom
    dvector N1(N);
    //Hamiltonian
    hermitian_matrix< complex<double>, upper> H(2 * L, 2 * L);
    //The matrix of P and Q coefficients of the fraction that represents green`s function`s matrix element on i-th atom
    matrix<double> greenFraction;
    //Number of iterations used
    unsigned int itr = 0;
    do
    { 
        ++itr;
        //Copy new values of the number d-electrons and magnetic moments and proceed the process 
        N[i] = N1[i];
        M[i] = M1[i];

        H.clear(); greenFraction.clear();
     
        formHamiltonian(H, tAngle, pAngle, N, M, E0, U0, hopingIntegrals, magneticField); 
               
        std::pair< dvector, dvector > cofficients = ThreeDiag(H, i, i, complex<double>(1, 0), complex<double>(1, 0));

        int ii0 = 0;
        greenFraction = Green(cofficients.first, cofficients.second, ii0);

        double Nu = CN(greenFraction, ii0);
        double Eu = CalculateEnergy(greenFraction, cofficients.first.size() - 1);      

        
        cofficients = ThreeDiag(H, i + L, i + L, complex<double>(1, 0), complex<double>(1, 0));
        greenFraction = Green(cofficients.first, cofficients.second, ii0);

        double Nd = CN(greenFraction, ii0);
        double Ed = CalculateEnergy(greenFraction, cofficients.first.size() - 1);

        cofficients = ThreeDiag(H, i, i + L, complex<double>(1, 0), complex<double>(1, 0)); 
        greenFraction = Green(cofficients.first, cofficients.second, ii0);
        double SFp = CN(greenFraction, ii0);

        cofficients = ThreeDiag(H, i, i + L, complex<double>(1, 0), complex<double>(-1, 0));
        greenFraction = Green(cofficients.first, cofficients.second, ii0);
        double SFn = CN(greenFraction, ii0);

        cofficients = ThreeDiag(H, i, i + L, complex<double>(0, 1), complex<double>(1, 0));
        greenFraction = Green(cofficients.first, cofficients.second, ii0);
        double AFp = CN(greenFraction, ii0);

        cofficients = ThreeDiag(H, i, i + L, complex<double>(0, 1), complex<double>(-1, 0));
        greenFraction = Green(cofficients.first, cofficients.second, ii0);
        double AFn = CN(greenFraction, ii0);

        N1[i] = Nu + Nd;
        M1[i] = (Nu - Nd) * cos(tAngle[i]) - ((SFp - SFn) * cos(pAngle[i]) - (AFp - AFn) * sin(pAngle[i])) * sin(tAngle[i]);

        E[i] = Eu + Ed - U0[i] * ( N1[i] * N1[i] - M1[i] * M1[i] ) * 0.25 ;
    }
    while ( (std::abs(M1[i] - M[i]) > delta || std::abs(N1[i] - N[i]) > delta)  && itr < maxIterCount ) ;
    N[i] = N1[i];
    M[i] = M1[i];
    
//    cout << i << ") N=" << 5 * N[i] << " M=" << 5 * M[i] << std::endl;
    if( itr == maxIterCount )
        throw "#Infinite selfconsist procedure.";
    return itr == 1;
}

void SConsistThreaded(int i, const std::pair< dvector, dvector >& angles, dvector& M, dvector& N, const dvector& E0, const dvector& U0, const dmatrix& hopingIntegrals, dvector& E, bool& isConsist)
{
    double magneticField = 0.0;
    //Saving magnetic moments for each atom
    dvector M1(M);
    //Saving number of d-electrons on each atom
    dvector N1(N);
    //Hamiltonian
//    hermitian_matrix< complex<double>, upper> H(2 * L, 2 * L);
    //The matrix of P and Q coefficients of the fraction that represents green`s function`s matrix element on i-th atom
//    matrix<double> greenFraction;
    //Number of iterations used
    unsigned int itr = 0;

    do
    {
        ++itr;
        //Copy new values of the number d-electrons and magnetic moments and proceed the process 
        N[i] = N1[i];
        M[i] = M1[i];

        hermitian_matrix< complex<double>, upper> H(2 * L, 2 * L);
        matrix<double> greenFraction;
        
        H.clear();
        greenFraction.clear();

        formHamiltonian(H, angles.first, angles.second, N, M, E0, U0, hopingIntegrals, magneticField); 

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
    while ( (std::abs(M1[i] - M[i]) > delta || std::abs(N1[i] - N[i]) > delta) && itr < maxIterCount ) ;
    N[i] = N1[i];
    M[i] = M1[i];
    
    if( itr == maxIterCount )
    {
        interrupt = true;
        return;
    }
    isConsist &= (itr == 1);
}

void correctErrors( dmatrix& results )
{
        for( int i = 0; i < results.size(); ++i )
        {
            std::vector< int > badDots;
            for( int j = 0; j < results.size(); ++j )
            {
                //We marked bad dots with -1 so if it is one then we remember its position
                if( results[i][j] == -1)
                {
                    badDots.push_back(j);
                    std::cout << "#Bad dot results[" << i << "][" << j << "]" << std::endl;
                }
                else
                {
                    while( !badDots.empty() )
                    {
                        double correction = results[i][ badDots.back() + 1 ];
                        int numCorrections = 1;
                        //Take into account dot over the bad one. It's 100% correct so no -1 check
                        if( i - 1 >= 0 )
                        {
                            correction += results[i - 1][badDots.back()];
                            numCorrections++;
                        }
                        //Take into account dot under the bad one
                        if( i + 1 < results.size() && results[i + 1][badDots.back()] != -1)
                        {
                            correction += results[i + 1][badDots.back()];
                            numCorrections++;
                        }
                        if( badDots.back() - 1 >= 0 && results[i][badDots.back() - 1] != - 1 )
                        {
                            correction += results[i][badDots.back() - 1];
                            numCorrections++;
                        }
                        correction /= numCorrections;

                        std::cout << "#results[" << i << "][" << badDots.back() << "] corrected with " << correction << std::endl; 
                        results[i][ badDots.back() ] = correction;
                        badDots.pop_back();
                    }
                }
            }          
            //Means that bad dots where till the end. exp: 2 3 5 6 -1 -1 -1 -1 (where -1 denotes a bad dot))
            if( !badDots.empty() )
            {
                // int firstGood = badDots[0] - 1;
                for (int j = badDots[0]; j < results.size(); ++j)
                {
                    int numCorrections = 0;
                    double correction = 0;
                    //Left to the bad dot
                    if( j - 1 >= 0 && results[i][j - 1] != -1)
                    {
                        correction += results[i][j - 1];
                        numCorrections++;
                    }
                    //Over the bad dot
                    if( i - 1 >= 0)
                    {
                        correction += results[i - 1][j];
                        numCorrections++;
                    }
                    //Under the bad dot
                    if( i + 1 < results.size() && results[i + 1][j] != -1)
                    {
                        correction += results[i + 1][j];
                        numCorrections++;
                    }                    
                    // results[i][j] = results[i][firstGood];
                    results[i][j] = correction / numCorrections;
                }
            }
        }    
}

dmatrix buildEnergySurface(dvector& thetaAngles, const dvector& phiAngles, dvector& magneticMoments, dvector& electronsNumber, 
                           const dvector& E0, const dvector& U0, const dmatrix& hopingIntegrals, const bool& threaded )
{
    std::vector< std::vector<double> > results;
    double step             = 0.1;
    double theta2Begin      = 0; //0.6;
    double theta2End        = 2 * 3.14;
    double theta3Begin      = 0;//0.6000;
    double theta3End        = 2 * 3.14; //0.6000;

    for( double theta2 = theta2Begin; theta2 <= theta2End; theta2 += step )
    {
        std::vector<double> bufResults;
        for( double theta3 = theta3Begin; theta3 <= theta3End;  theta3 += step )
        {
            try
            {
                //Energy
                std::vector<double> Energy(L, 0); 

                std::vector<double> N(electronsNumber);
                //initial magnenic moments
                std::vector<double> M(magneticMoments);

                thetaAngles[0] = 0;
                thetaAngles[1] = theta2;
                thetaAngles[2] = theta3;

                bool isConsist = true;
                unsigned int iterations = 0;

                std::cout << "#" << theta2 << " " << theta3 << std::endl;

                if( !threaded )
                {     
                    do
                    {
                        isConsist = true;
                        for (int i = 0; i < L; ++i)
                            isConsist &= SConsist(i, thetaAngles, phiAngles, M, N, E0, U0, hopingIntegrals, Energy);
                        if( ++iterations == 50 )
                                throw "#Infinite consist cluster loop.";    
                    }
                    while (!isConsist);

                }
                else
                {
                    do
                    {
                            int j = 0;
                            isConsist = true;
                            while( j < L )
                            {	  
                                    boost::thread_group threads;
                                    for ( int i = 0; i < CORE_COUNT && j < L; ++i, ++j )
                                            threads.add_thread( new boost::thread( &SConsistThreaded, j,  std::make_pair(thetaAngles, phiAngles), boost::ref(M), boost::ref(N), E0, U0, hopingIntegrals, boost::ref(Energy), boost::ref(isConsist)  ) );
                                    threads.join_all();
                                    if( interrupt == true )
                                    {
                                        interrupt = false;
                                        throw "#Infinite selfconsist procedure.";  
                                    }
                            }
                            if( ++iterations == 10 )
                                throw "#Infinite consist cluster loop.";    

                    } while (!isConsist);
                }
//                GradE(thetaAngles, phiAngles, M, N, E0, U0, hopingIntegrals);

                bufResults.push_back( std::accumulate( Energy.begin(), Energy.end(), 0.0 ) );

                std::cout << "#Number of d-electrons: ";
                for (int i = 0; i < N.size(); ++i)
                    std::cout << 5 * N[i] << " ";
                std::cout << std::endl;                    
/*
                std::cout << "#Magnetic monents: ";
                for (int i = 0; i < M.size(); ++i)
                    std::cout << 5 * M[i] << " ";
                std::cout << std::endl;  
 *                   
*/
            }
            catch (const char* msg)
            {
                std::cout << msg;
//                 std::cout << "#No solution: " << msg <<  std::endl;
                 std::terminate();
                bufResults.push_back(-1);
                continue;
            }
            catch( boost::thread_interrupted const& e )
            {
//                 std::cout << "#No solution" <<  std::endl;
                std::terminate();
                bufResults.push_back(-1);
                continue;                    
            }
        }                

        results.push_back( bufResults );
    }

    //Approximates the point where selfconsistent result was not found
//    correctErrors(results);

    return results;
}

dvector GradE(dvector& thetaAngles, const dvector& phiAngles, dvector& magneticMoments, dvector& electronsNumber, 
                           const dvector& E0, const dvector& U0, const dmatrix& hopingIntegrals)
{
    cx_mat H(2 * L, 2 * L);
    H.fill(0);
    formHamiltonian(H, thetaAngles, phiAngles, electronsNumber, magneticMoments, E0, U0, hopingIntegrals, 0);    
    //Make the matrix Hermitian
//    H += H.t();
   
    //Finding the basis in which H is a diagonal transformation by calculating eigenvalues
    cx_mat eigenVectors;
    vec eigenValues;
    eig_sym(eigenValues, eigenVectors, H);
    
    //2 * L means we have L atoms with two degrees of freedom. First L terms terms 
    dvector gradient(2 * L);
    
   
    for(unsigned int i = 0; i < gradient.size(); ++i)
    {
        //H derivative with respect to theta angles
        cx_mat dHdt(2 * L, 2 * L);
        //H derivative with respect to phi angles        
        cx_mat dHdp(2 * L, 2 * L);
        
        dHdt.fill(0);
        dHdp.fill(0);

        if(i < L)
        {
            dHdt(i, i)          = complex<double>( U0[i] * 0.5 * magneticMoments[i] * sin(thetaAngles[i]), 0);
            dHdt(i + L, i + L)  = complex<double>(-U0[i] * 0.5 * magneticMoments[i] * sin(thetaAngles[i]), 0); 
            double x            = U0[i] * magneticMoments[i] * cos(thetaAngles[i]) * 0.5;
            dHdt(i, i + L)      = complex<double>(x * cos(phiAngles[i]), -x * sin(phiAngles[i]) );
            dHdt(i + L, i)      = complex<double>(x * cos(phiAngles[i]),  x * sin(phiAngles[i]) );
            
            //Moving dHdt to the basis of our Hamiltonian. H**T equals H inverse in case of such basis
            cx_mat dHdt_newBasis = eigenVectors.t() * dHdt * eigenVectors;
                                  
            double s1 = 0;
            double s2 = 0;

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
            
            double s1 = 0;
            double s2 = 0;

            for(unsigned p = 0; p < 2 * L; ++p)
            {
                
                s1 += dHdp_newBasis(p, p).real();
                s2 += dHdp_newBasis(p, p).real() * atan( eigenValues(p) );
            }
        }
        gradient[i] = (0.5 * s1 - s2 / boost::math::constants::pi<double>() );
       
    }
    
    return gradient;
}

void buildSelfconsistentSolution(const int& imageNum, const mat& angles, dmatrix& magneticMoments, dmatrix& electronsNumber, 
                           const dvector& E0, const dvector& U0, const dmatrix& hopingIntegrals, dmatrix& Gradients, dvector& Energies)
{
    bool isConsist = false;
    unsigned int iterations = 0;

    dvector thetaAngles;
    dvector phiAngles;
    dvector Energy(L, 0);
    dvector N(electronsNumber[imageNum]);
    dvector M(magneticMoments[imageNum]); 
    for(int i = 0; i < L; ++i)
    {
        thetaAngles.push_back( angles(imageNum, i));
        phiAngles.push_back( angles(imageNum, i + L));
    }
    
    do
    {
        isConsist = true;
        for (int i = 0; i < L; ++i)
            isConsist &= SConsist(i, thetaAngles, phiAngles, M, N, E0, U0, hopingIntegrals, Energy);
        if( ++iterations == 50 )
                throw "#Infinite consist cluster loop.";    
    }
    while (!isConsist);
     
    electronsNumber[imageNum]   = N;
    magneticMoments[imageNum]   = M;
    Gradients[imageNum]         = GradE(thetaAngles, phiAngles, M, N, E0, U0, hopingIntegrals);
    Energies[imageNum]          =  std::accumulate( Energy.begin(), Energy.end(), 0.0 );
}

dvector pathTangent(const int& imageNum, const mat& path, const dvector& energies)
{
    double E1 = energies[imageNum - 1];
    double E2 = energies[imageNum];
    double E3 = energies[imageNum + 1];
    dvector tau(2 * L);
    
//    double length = 0;
//    for(int i = 0; i < 2 * L; ++i)
//    {
//        tau[i] = (path(imageNum + 1, i) - path(imageNum - 1, i) );
//        length += tau[i] * tau[i];
//    }
//    for(int i = 0 ; i < 2 * L; ++i)
//    {
//        tau[i] /= std::sqrt(length);
//    }
//        return tau;
    
    
    
    /*
    if( E3 > E1 )
    {
        for(int i = 0; i < 2 * L; ++i)
            tau[i] = std::max( std::abs(E3 - E2), std::abs(E2 - E1)) * ( path(imageNum + 1, i) - path(imageNum, i) ) + std::min( std::abs(E3 - E2), std::abs(E2 - E1)) * ( path(imageNum , i) - path(imageNum - 1, i) );
    }    
    else if( E3 < E1 )
        for(int i = 0; i < 2 * L; ++i)
            tau[i] = std::min( std::abs(E3 - E2), std::abs(E2 - E1)) * ( path(imageNum + 1, i) - path(imageNum, i) ) + std::max( std::abs(E3 - E2), std::abs(E2 - E1)) * ( path(imageNum , i) - path(imageNum - 1, i) );
    else
        for(int i = 0; i < 2 * L; ++i)
            tau[i] = (path(imageNum + 1, i) - path(imageNum - 1, i) );
     */
    if( E3 > E2 && E2 > E1)
        for(int i = 0; i < 2 * L; ++i)
            tau[i] = path(imageNum + 1, i) - path(imageNum, i);
    else if( E3 < E2 && E2 < E1)
        for(int i = 0; i < 2 * L; ++i)
            tau[i] = path(imageNum, i) - path(imageNum - 1, i);
    else
    {
        double dEmax = std::abs(E3 - E2);
        double dEmin = std::abs(E1 - E2);

        if(dEmax < dEmin)
            std::swap(dEmax, dEmin);
        if( E3 > E1 )
                for(int i = 0; i < 2 * L; ++i)
                    tau[i] = dEmax * ( path(imageNum + 1, i) - path(imageNum, i) ) + dEmin * ( path(imageNum, i) - path(imageNum - 1, i) );
        else
                for(int i = 0; i < 2 * L; ++i)
                    tau[i] = dEmin * ( path(imageNum + 1, i) - path(imageNum, i) ) + dEmax * ( path(imageNum, i) - path(imageNum - 1, i) );            
    }

    double length = 0;
        for(int i = 0 ; i < 2 * L; ++i)
            length += tau[i] * tau[i];
    for(int i = 0 ; i < 2 * L; ++i)
        tau[i] /= std::sqrt(length);
    return tau;

}
void buldMEP(dvector& thetaAngles, const dvector& phiAngles, dvector& magneticMoments, dvector& electronsNumber, 
                           const dvector& E0, const dvector& U0, const dmatrix& hopingIntegrals, int numberOfImages, const dvector& initialState, const dvector& finalState )
{
    dmatrix pathN(numberOfImages);
    dmatrix pathM(numberOfImages);    
    dmatrix gradients(numberOfImages);
    dvector energies(numberOfImages);
    mat path(numberOfImages, 2 * L);

    for(int i = 0; i < 2 * L; ++i)
    {
        double dh = (finalState[i] - initialState[i]) / (numberOfImages - 1);
        for(int j = 0; j < numberOfImages; ++j)
            path(j, i) = initialState[i] + dh * j;
    }  
    
    for(int i = 0; i < numberOfImages; ++i)
    {
        for(int j = 0; j < L; ++j)
        {
            pathN[i].push_back( electronsNumber[j] );
            pathM[i].push_back( magneticMoments[j] );
        }
    }
    
    buildSelfconsistentSolution(0                 , path, pathM, pathN, E0, U0,  hopingIntegrals, gradients, energies);
    buildSelfconsistentSolution(numberOfImages - 1, path, pathM, pathN, E0, U0,  hopingIntegrals, gradients, energies);

    double forceTol = 0.1;

    for(int i = 0; i < numberOfImages; ++i )
    {
        cout << path(i, 1) * 60 / 2 / 3.14 << " " << path(i, 2) * 60 / 2 / 3.14  << " " << 0 << std::endl;
    }
    do
    {
        for(int i = 1; i < numberOfImages -1; ++i)
                buildSelfconsistentSolution(i, path, pathM, pathN, E0, U0,  hopingIntegrals, gradients, energies); 

        double totalForce = 0;
        for(int i = 1; i < numberOfImages - 1; ++i)
        {
            dvector tauVector = pathTangent(i, path, energies);
            
            dvector springForce(2 * L, 0);
            for(int j = 1; j < 2 * L; ++j)
                springForce[j]     =  0.5 * ( path(i + 1, j) + path(i - 1, j) - 2 * path(i, j) );
            double gradientProjection    = 0;
            double springForceProjection = 0;
            for( int j = 1; j < 2 * L; ++j )
            {
                gradientProjection      += gradients[i][j] * tauVector[j];
                springForceProjection   += springForce[j] * tauVector[j];
            }
            
            for( int j = 1; j < 2 * L; ++j )
            {
                path(i, j) += -( gradients[i][j] - tauVector[j] * gradientProjection) + springForceProjection * tauVector[j]; 
                totalForce +=  ( gradients[i][j] - tauVector[j] * gradientProjection) * ( gradients[i][j] - tauVector[j] * gradientProjection);
            }
        }
        cout << "Total force: " << totalForce << std::endl;
        if( std::abs(totalForce)  < 1e-5 )
            break;
    } while( true );

    cout << "Final path:\n";
    for(int i = 0; i < numberOfImages; ++i )
    {
        cout << path(i, 1) * 60 / 2 / 3.14 << " " << path(i, 2) * 60 / 2 / 3.14  << " " << 0 << std::endl;
    }
    
}

int main(int argc, char* argv[])
{
    std::cout << std::setprecision(16);
        //Hoping itegrals between atoms
    dmatrix hoping_integrals;
    //Position of non-perturbed d-level with respect to the Fermi energy for all atoms
    std::vector<double> E0;
    //Coulomb repulsion energy for all atoms
    std::vector<double> U0;
    //initial number of d-electrons
    std::vector<double> ElectronsNumber;
    //initial magnenic moments
    std::vector<double> MagneticMoments;
    //Theta angles for atoms
    std::vector<double> thetaAngles;
    //Phi angle for atoms
    std::vector<double> phiAngles;
    //Magnetic field
    double dh = 0;
    
    try 
    {
        std::ifstream hparams_file("hoping.txt");
        std::string buf;
        while (std::getline(hparams_file, buf))
        {
            std::istringstream line(buf);
            hoping_integrals.push_back(std::vector<double>());
            std::copy(std::istream_iterator<double>(line), std::istream_iterator<double>(), std::back_inserter(hoping_integrals.back()));
        }
        //Reading all necesary params from the text file
        std::ifstream init_file("InitFile.txt");

        GetParam(init_file, E0);
        GetParam(init_file, U0);
        GetParam(init_file, ElectronsNumber);
        GetParam(init_file, MagneticMoments);
        GetParam(init_file, thetaAngles);
        GetParam(init_file, phiAngles);
            
//        TestGradient(thetaAngles, phiAngles, MagneticMoments, ElectronsNumber, E0, U0, hoping_integrals);
        buldMEP(thetaAngles, phiAngles, MagneticMoments, ElectronsNumber, E0, U0, hoping_integrals, 8, {0, 0.2 * 3.14, 1.2 * 3.14, 0, 0, 0}, {0, 1.5 * 3.14, 0.5 * 3.14, 0, 0, 0});
        return 0;
        dmatrix results = buildEnergySurface(thetaAngles, phiAngles, MagneticMoments, ElectronsNumber, E0, U0, hoping_integrals, false );
 
        std::ofstream logFile;
        logFile.open( "temp.txt", std::ofstream::out | std::ofstream::trunc );
        for( int i = 0; i < results.size(); ++i )
        {
            for( int j = 0; j < results.size(); ++j )
                logFile << results[i][j] << " ";
            logFile <<std::endl;
        }
    }
    
            
    catch (std::exception& exp)
    {
        std::cout << exp.what();
    }
    catch( const char* msg)
    {
        cout << msg;
    }
    return 0;
}
