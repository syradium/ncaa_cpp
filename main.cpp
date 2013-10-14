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

const int CORE_COUNT = 4;

bool interrupt = false;

const int L = 3;
//The maximum number of iteration in selfconsist procedure. If that number is reached we consider that selfconsist procedure will never end, so we terminate it.
const int maxIterCount = 100;

using std::complex;
using std::numeric_limits;
using namespace boost::numeric::ublas;
using namespace boost::math::tools;

typedef std::vector<std::vector<double> > dmatrix;
typedef std::vector< double > dvector;
typedef hermitian_matrix<std::complex<double>, upper> herm_matrix;

const double delta = 1e-8 + std::numeric_limits<double>::epsilon();

//H - hamiltonian, tAngle - thetta angles, pAngle - phi angles, N - number of d-electrons on each atom, M - magnetic moments of each atom
void formHamiltonian(herm_matrix& H, const dvector& tAngle, const dvector& pAngle, const dvector& N, const dvector& M, const dvector& E0, const dvector& U0, const dmatrix& V, const double& dh)
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
    for (imax = 0; imax < (2 * L - 1) && x * x > 1e-8; ++imax)
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

double f(const double& x, const double& a0, const double& b0, const std::vector<double>& p0, const std::vector<double>& q0)
{
    double s = 0;
    for (int i = 0; i < p0.size() /*&& std::abs(p0[i]) >= 1E-15*/ /*0.00000001*/; ++i)
        if( std::abs(p0[i]) >= 1E-15 )
                s += p0[i] / (x - q0[i]);
    return x - a0 - b0 * s;
}

double bisection(const double& a, const double& b, const double& a0, const double& b0, const std::vector<double>& p0, const std::vector<double>& q0)
{

    double l = a;
    double r = b;
    const double EPS = delta * 10;
    while (r - l > EPS)
    {
        double m = (l + r) / 2;
        
        if( f(m, a0, b0, p0, q0) * f(r, a0, b0, p0, q0) <= 0 )
            l = m;
        else
            r = m;
//        if (f(m, a0, b0, p0, q0) > 0 && f(r, a0, b0, p0, q0) > 0 || f(m, a0, b0, p0, q0) < 0 && f(r, a0, b0, p0, q0) < 0)
//            r = m;
//        else
//            l = m;
    }
    return l;

}

double bisection2(const double& a, const double& b, const double& a0, const double& b0, const std::vector<double>& p0, const std::vector<double>& q0)
{
    double l = a;
    double r = b;
    const double EPS = 1E-10;
    unsigned iter = 0;
    while (r - l > EPS)
    {
        double m = (l + r) / 2;

//        std::cout << m << " " << f(m, a0, b0, p0, q0) << std::endl;
        if(++iter > 50)
            std::terminate();
        if( f(m, a0, b0, p0, q0) * f(r, a0, b0, p0, q0) <= 0 )
            l = m;
        else
            r = m;
//        if (f(m, a0, b0, p0, q0) > 0 && f(r, a0, b0, p0, q0) > 0 || f(m, a0, b0, p0, q0) < 0 && f(r, a0, b0, p0, q0) < 0)
//            r = m;
//        else
//            l = m;
    }
    return l;

}




void FR(const double& a, const double& b, matrix<double>& t, std::vector<double>& p0, std::vector<double>& q0)
{
    double a0 = a;
    double b0 = b; //* b;
    std::vector<double> z(t.size2(), 0);
    int ii0 = 0;

    for (int i = 2 * L - 1; i >= 0; --i)
    {
        p0[i] = t(0, i); 
        q0[i] = t(1, i);
        if (std::abs(p0[i] * b0) <= 1e-10)
            ii0 = i - 1;
    }

    if (ii0 < 0)
    {
        t(0, 0) = 1;
        t(1, 0) = a;
        return;
    }
   
//    double leftBound    = q0[0] - 10 - std::abs(a);
//    double rightBound   = q0[ii0] + 1 + std::abs(a);
//    for( int i = 0; i < p0.size() ; ++i )
//    {
//        leftBound -=  p0[i];
//        rightBound += p0[i];
//    }
    double leftBound = q0[0]- 1E4;
    double rightBound = q0[ii0] + 1E4;
    
    z[0] = bisection( leftBound, q0[0] - delta, a0, b0, p0, q0);

    for (int i = 0; i <= ii0; ++i)
        z[i + 1] = bisection(q0[i] + delta, q0[i + 1] - delta, a0, b0, p0, q0);
    z[ii0 + 1] = bisection(q0[ii0] + delta, rightBound, a0, b0, p0, q0);
    
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
}


//Calculates the coefficients P and Q of the Green`s function`s matrix element represented by a fraction
//"a" and "b" are the diagonal and subdiagonal elemnts of the hamiltonian

matrix<double> Green(const std::vector<double>& a, const std::vector<double>& b)
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
    
    for (int i = imax - 1; i >= 0; --i)
        FR(a[i], b[i], t, p0, q0);
    return t;
}

void FR2(const double& a, const double& b, matrix<double>& t, std::vector<double>& p0, std::vector<double>& q0)
{
    double a0 = a;
    double b0 = b; //* b;
    std::vector<double> z(t.size2(), 0);
    int ii0 = 0;

    for (int i = 2 * L - 1; i >= 0; --i)
    {
        p0[i] = t(0, i); 
        q0[i] = t(1, i);
        if (std::abs(p0[i] * b0) <= 1E-10)//<= 0.00000001)        
            ii0 = i - 1;
    }

    if (ii0 < 0)
    {
        t(0, 0) = 1;
        t(1, 0) = a;
        return;
    }
   
    double leftBound    = q0[0] - 10 - std::abs(a);
    double rightBound   = q0[ii0] + 1 + std::abs(a);
    for( int i = 0; i < p0.size() ; ++i )
    {
        leftBound -=  p0[i];
        rightBound += p0[i];
    }
//    double leftBound = q0[0]- 1E4;
//    double rightBound = q0[ii0] + 1E4;
    
    
    z[0] = bisection2( leftBound, q0[0] - 1e-10, a0, b0, p0, q0);
    std::cout << "leftBound=" << leftBound << std::endl;
    std::cout << "z[0]" << z[0] << std::endl;

    for (int i = 0; i <= ii0; ++i)
        z[i + 1] = bisection(q0[i] + 1e-10, q0[i + 1] - 1e-10, a0, b0, p0, q0);
    z[ii0 + 1] = bisection(q0[ii0] + 1e-11, rightBound, a0, b0, p0, q0);
    
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
}



matrix<double> Green2(const std::vector<double>& a, const std::vector<double>& b)
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
    
    for (int i = imax - 1; i >= 0; --i)
        FR2(a[i], b[i], t, p0, q0);
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

double CN2(const matrix<double>& t, const int& imax)
{
    double s1 = 0;
    double s2 = 0;
    for (int i = 0; i <= imax; ++i)
    {
        s1 += t(0, i);
        s2 += t(0, i) * std::atan(t(1, i));
    }
    std::cout << "sum(s1): " << s1 << std::endl;
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
               
        greenFraction = Green(cofficients.first, cofficients.second);

        double Nu = CN(greenFraction, cofficients.first.size() - 1);
        double Eu = CalculateEnergy(greenFraction, cofficients.first.size() - 1);      

        
        cofficients = ThreeDiag(H, i + L, i + L, complex<double>(1, 0), complex<double>(1, 0));
        greenFraction = Green(cofficients.first, cofficients.second);

        double Nd = CN(greenFraction, cofficients.first.size() - 1);
        double Ed = CalculateEnergy(greenFraction, cofficients.first.size() - 1);
      
        cofficients = ThreeDiag(H, i, i + L, complex<double>(1, 0), complex<double>(1, 0));
 
        greenFraction = Green(cofficients.first, cofficients.second);

        double SFp = CN(greenFraction, cofficients.first.size() - 1);

        cofficients = ThreeDiag(H, i, i + L, complex<double>(1, 0), complex<double>(-1, 0));
        greenFraction = Green(cofficients.first, cofficients.second);

        double SFn = CN(greenFraction, cofficients.first.size() - 1);

        cofficients = ThreeDiag(H, i, i + L, complex<double>(0, 1), complex<double>(1, 0));
        greenFraction = Green(cofficients.first, cofficients.second);
        double AFp = CN(greenFraction, cofficients.first.size() - 1);

        cofficients = ThreeDiag(H, i, i + L, complex<double>(0, 1), complex<double>(-1, 0));
        greenFraction = Green(cofficients.first, cofficients.second);
        double AFn = CN(greenFraction, cofficients.first.size() - 1);
       
        N1[i] = Nu + Nd;
        M1[i] = (Nu - Nd) * cos(tAngle[i]) - ((SFp - SFn) * cos(pAngle[i]) - (AFp - AFn) * sin(pAngle[i])) * sin(tAngle[i]);
        E[i] = Eu + Ed - U0[i] * ( N1[i] * N1[i] - M1[i] * M1[i] ) * 0.25 ;
    }
    while ( (std::abs(M1[i] - M[i]) > delta || std::abs(N1[i] - N[i]) > delta) && itr < maxIterCount ) ;
    N[i] = N1[i];
    M[i] = M1[i];
    if( itr == maxIterCount )
        throw "#Infinite selfconsist procedure.";
    return itr == 1;
}

void SConsistThreaded(int i, const std::pair< dvector, dvector >& angles, /*const dvector& tAngle, const dvector& pAngle, */dvector& M, dvector& N, const dvector& E0, const dvector& U0, const dmatrix& hopingIntegrals, dvector& E, bool& isConsist)
{
    double magneticField = 0.08;
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
        
        
        greenFraction = Green(cofficients.first, cofficients.second);
        double Nu = CN(greenFraction, cofficients.first.size() - 1);
        double Eu = CalculateEnergy(greenFraction, cofficients.first.size() - 1);      

        
        cofficients = ThreeDiag(H, i + L, i + L, complex<double>(1, 0), complex<double>(1, 0));
        greenFraction = Green(cofficients.first, cofficients.second);

        double Nd = CN(greenFraction, cofficients.first.size() - 1);
        double Ed = CalculateEnergy(greenFraction, cofficients.first.size() - 1);
      
        cofficients = ThreeDiag(H, i, i + L, complex<double>(1, 0), complex<double>(1, 0));
 
        greenFraction = Green(cofficients.first, cofficients.second);

        double SFp = CN(greenFraction, cofficients.first.size() - 1);

        cofficients = ThreeDiag(H, i, i + L, complex<double>(1, 0), complex<double>(-1, 0));
        greenFraction = Green(cofficients.first, cofficients.second);

        double SFn = CN(greenFraction, cofficients.first.size() - 1);

        cofficients = ThreeDiag(H, i, i + L, complex<double>(0, 1), complex<double>(1, 0));
        greenFraction = Green(cofficients.first, cofficients.second);
        double AFp = CN(greenFraction, cofficients.first.size() - 1);

        cofficients = ThreeDiag(H, i, i + L, complex<double>(0, 1), complex<double>(-1, 0));
        greenFraction = Green(cofficients.first, cofficients.second);
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
        // std::cout << "#Interruption" << std::endl;
        return;
    }
    isConsist &= (itr == 1);
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

                    std::vector<double> N(ElectronsNumber);
                    //initial magnenic moments
                    std::vector<double> M(MagneticMoments);

                    thetaAngles[0] = 0;
                    thetaAngles[1] = theta2;
                    thetaAngles[2] = theta3;

                    bool isConsist = true;
                    unsigned int iterations = 0;

                    std::cout << "#" << theta2 << " " << theta3 << std::endl;

  /*                  
                    do
                    {
                        isConsist = true;
                        for (int i = 0; i < L; ++i)
                            isConsist &= SConsist(i, thetaAngles, phiAngles, M, N, E0, U0, hoping_integrals, Energy);
                        if( ++iterations == 10 )
                                throw "#Infinite consist cluster loop.";    
                    }
                    while (!isConsist);
*/                  
                    do
                    {
                            int j = 0;
                            isConsist = true;
                            while( j < L )
                            {	  
                                    boost::thread_group threads;
                                    for ( int i = 0; i < CORE_COUNT && j < L; ++i, ++j )
                                            threads.add_thread( new boost::thread( &SConsistThreaded, j,  std::make_pair(thetaAngles, phiAngles), boost::ref(M), boost::ref(N), E0, U0, hoping_integrals, boost::ref(Energy), boost::ref(isConsist)  ) );
                                    threads.join_all();
                                    if( interrupt == true )
                                    {
                                        interrupt = false;
                                        throw "#Infinite selfconsistconsist procedure.";  
                                    }
                            }
                            if( ++iterations == 10 )
                                throw "#Infinite consist cluster loop.";    

                    } while (!isConsist);

                    double energy = 0;
                    for( int i = 0; i < Energy.size(); ++i)
                        energy += Energy[i];
                    bufResults.push_back( energy );
                    
/*                    std::cout << "#Number of d-electrons: ";
                    for (int i = 0; i < N.size(); ++i)
                        std::cout << 5 * N[i] << " ";
                    std::cout << std::endl;                    

                    std::cout << "#Magnetic monents: ";
                    for (int i = 0; i < M.size(); ++i)
                        std::cout << 5 * M[i] << " ";
                    std::cout << std::endl;                    
*/
                }
                catch (const char* msg)
                {
                    // std::cout << "#No solution: " << msg <<  std::endl;
                    bufResults.push_back(-1);
                    continue;
                }
                catch( boost::thread_interrupted const& e )
                {
                    // std::cout << "#No solution" <<  std::endl;
                    bufResults.push_back(-1);
                    continue;                    
                }
                
                /*		do
                               {
                                       int j = 0;
                                       isConsist = true;
                                       while( j < L )
                                       {	  
                                               boost::thread_group threads;
                                               for ( int i = 0; i < CORE_COUNT && j < L; ++i, ++j )
                                               {	
                                                       threads.add_thread( new boost::thread( &SConsist, j, t, p, boost::ref(M), boost::ref(N), E0, U0, hoping_integrals, boost::ref(isConsist)  ) );
                                                       std::cout << "Thread launched" << std::endl;	
                                               }
                                               threads.join_all();
                                               std::cout << "Block end" << std::endl;
                                       }
                               } while (!isConsist);
                */
            }
            results.push_back( bufResults );
        }
        
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
            //Means that bad dots where till the end. exp: 2 3 5 6 -1 -1 -1 -1
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

        
        for( int i = 0; i < results.size(); ++i )
        {
            for( int j = 0; j < results.size(); ++j )
                std::cout << results[i][j] << " ";
            std::cout << std::endl;
        }
    }
    catch (std::exception& exp)
    {
        std::cout << exp.what();
    }
    return 0;
}
