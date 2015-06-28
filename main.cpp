/* 
 * File:   main.cpp
 * Author: Radium
 *
 * Created on January 14, 2013, 1:15 AM
 */

#define BOOST_LOG_DYN_LINK 1

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
//#include <c++/4.7/complex>
//#include <c++/4.7/bits/stl_numeric.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

const int CORE_COUNT = 1;

bool interrupt = false;

const int L = 3;
//The maximum number of iteration in selfconsist procedure. If that number is reached we consider that selfconsist procedure will never end, so we terminate it.
const int maxIterCount = 400;

using std::complex;
using std::numeric_limits;
using namespace boost::numeric::ublas;
using namespace boost::math::tools;
using namespace arma;
namespace logging = boost::log;


typedef std::vector<std::vector<double> > dmatrix;
typedef std::vector< double > dvector;
typedef hermitian_matrix<std::complex<double>, upper> herm_matrix;


const double delta  = 1e-10;
const double eps    = 1e-11 + std::numeric_limits<double>::epsilon(); 

//#include "debug_functions.h"


void recordPath( const mat& path, const char* fileName);
void recordPoint( const mat& path, const dmatrix& pathN, const dmatrix& pathM, const dvector& energies, const char* fileName, const int& pointIndex );
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

		std::cout << N1[i] <<  " " << M1[i] << std::endl; 
		if(itr == 5)
			throw "";
		E[i] = Eu + Ed - U0[i] * ( N1[i] * N1[i] - M1[i] * M1[i] ) * 0.25 ;
	}
	while ( (std::abs(M1[i] - M[i]) > delta || std::abs(N1[i] - N[i]) > delta)  && itr < maxIterCount ) ;
	N[i] = N1[i];
	M[i] = M1[i];

	if( itr == maxIterCount )
		throw "#Infinite selfconsist procedure.";
	return itr == 1;
}

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

dmatrix buildEnergySurface(dvector thetaAngles, dvector phiAngles, const dvector& magneticMoments, const dvector& electronsNumber, 
		const dvector& E0, const dvector& U0, const dmatrix& hopingIntegrals, const bool& threaded )
{
	std::vector< std::vector<double> > results;
	int stepNumber          =  99;
	double theta2Begin      = -1.0 * boost::math::constants::pi<double>();
	double theta2End        =  3.0 * boost::math::constants::pi<double>();
	double theta3Begin      = -0.0 * boost::math::constants::pi<double>();
	double theta3End        =  4.0 * boost::math::constants::pi<double>();

	BOOST_LOG_TRIVIAL(info) << "Building energy surface with " << stepNumber << " steps";
	for( int th2 = 0; th2 < stepNumber; ++th2 )
	{
		BOOST_LOG_TRIVIAL(info) << "Step th2 " << th2 + 1 << " of " << stepNumber;

		std::vector<double> bufResults;
		thetaAngles[0] = 0;
		thetaAngles[1] = theta2Begin + (theta2End - theta2Begin) * th2 / (stepNumber - 1);

		for( int th3 = 0; th3 < stepNumber; ++th3 )
		{
			try
			{
				//Energy
				std::vector<double> Energy(L, 0); 

				std::vector<double> N(electronsNumber);
				//initial magnenic moments
				std::vector<double> M(magneticMoments);

				thetaAngles[2] = theta3Begin + (theta3End - theta3Begin) * th3 / (stepNumber - 1);

				bool isConsist = true;
				unsigned int iterations = 0;

				BOOST_LOG_TRIVIAL(debug) << "#" << thetaAngles[1] << " " << thetaAngles[2];

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

							std::pair< dvector, dvector > resultNM = std::make_pair(N, M);
							for ( int i = 0; i < CORE_COUNT && j < L; ++i, ++j )
							{
								threads.add_thread( new boost::thread( &SConsistThreaded, j, boost::ref(resultNM), std::make_pair(thetaAngles, phiAngles), 
											M, N, std::make_pair(E0, U0), 
											hopingIntegrals, boost::ref(Energy), boost::ref(isConsist)  ) );
							}
							threads.join_all();
							if( interrupt == true )
							{
								interrupt = false;
								throw "#Infinite selfconsist procedure.";  
							}
							N = resultNM.first;
							M = resultNM.second;
						}
						//                            if( ++iterations == 500 )
						//                                throw "#Infinite consist cluster loop.";    

					} while (!isConsist);
				}
				//                GradE(thetaAngles, phiAngles, M, N, E0, U0, hopingIntegrals);

				bufResults.push_back( std::accumulate( Energy.begin(), Energy.end(), 0.0 ) );

				BOOST_LOG_TRIVIAL(debug) << "#Number of d-electrons: ";
				for (int i = 0; i < N.size(); ++i)
					BOOST_LOG_TRIVIAL(debug)  << 5 * N[i] << " ";

				/*/
				  std::cout << "#Magnetic monents: ";
				  for (int i = 0; i < M.size(); ++i)
				  std::cout << 5 * M[i] << " ";
				  std::cout << std::endl;                    
				  */
			}
			catch (const char* msg)
			{
				BOOST_LOG_TRIVIAL(error) << msg;
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

	double energyMinima = results[0][0];
	for(int i = 0; i < results.size(); ++i)
		for(int j = 0; j < results[i].size(); ++j)
			if( results[i][j] < energyMinima )
				energyMinima = results[i][j];
	for(int i = 0; i < results.size(); ++i)
		for(int j = 0; j < results[i].size(); ++j)
			results[i][j] -= energyMinima;

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

void buildSelfconsistentSolution(const int& imageNum, const mat& angles, dmatrix& magneticMoments, dmatrix& electronsNumber, 
		const dvector& E0, const dvector& U0, const dmatrix& hopingIntegrals, dmatrix& Gradients, dvector& Energies)
{
	bool isConsist = false;
	unsigned int iterations = 0;

	dvector thetaAngles(L, 0);
	dvector phiAngles(L, 0);
	dvector Energy(L, 0);
	dvector N(electronsNumber[imageNum]);
	dvector M(magneticMoments[imageNum]); 
	for(int i = 0; i < L; ++i)
	{
		thetaAngles[i] = angles(imageNum, i);
		phiAngles[i]   = angles(imageNum, i + L);
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
	Energies[imageNum]          = std::accumulate( Energy.begin(), Energy.end(), 0.0 );
}

void buildSelfCSolution(dvector& thetaAngles, const dvector& phiAngles, dvector& magneticMoments, dvector& electronsNumber, 
		const dvector& E0, const dvector& U0, const dmatrix& hopingIntegrals, dvector& Energy)
{
	bool isConsist = false;
	unsigned int iterations = 0;

	do
	{
		isConsist = true;
		for (int i = 0; i < L; ++i)
			isConsist &= SConsist(i, thetaAngles, phiAngles, magneticMoments, electronsNumber, E0, U0, hopingIntegrals, Energy);
		if( ++iterations == 50 )
			throw "#Infinite consist cluster loop.";    
	}
	while (!isConsist);
}

//dmatrix hessian(dvector& thetaAngles, const dvector& phiAngles, dvector& magneticMoments, dvector& electronsNumber, 
//                           const dvector& E0, const dvector& U0, const dmatrix& hopingIntegrals)
//{
//    dvector Energy(L, 0); 
//    dvector N(electronsNumber);
//    dvector M(magneticMoments);
//    cx_mat H(2 * L, 2 * L);
//    double dx = 0.4;
//    buildSelfCSolution( thetaAngles, phiAngles, M, N, E0, U0, hopingIntegrals, Energy);
//
//    double en = std::accumulate( Energy.begin(), Energy.end(), 0.0 );
//
//    for(int i = 0; i < 2 * L; ++i)
//    {
//        double tp1 = thetaAngles[i];
//        
//        Energy.clear();
//        N = electronsNumber;
//        M = magneticMoments;
//        
//        thetaAngles[i] = tp1 + dx;
//        buildSelfCSolution( thetaAngles, phiAngles, M, N, E0, U0, hopingIntegrals, Energy);
//        double en1 = std::accumulate( Energy.begin(), Energy.end(), 0.0 );
//
//        thetaAngles[i] = tp1 - dx;
//        buildSelfCSolution( thetaAngles, phiAngles, M, N, E0, U0, hopingIntegrals, Energy);
//        double en2 = std::accumulate( Energy.begin(), Energy.end(), 0.0 );
//        
//        H(i, i) = (en1 + en2 - 2 * en) / (dx * dx);
//        for(int j = i + 1; j < 2 * L; ++j)
//        {
//                        
//        }
//    }
//}

dvector pathTangent(const int& imageNum, const mat& path, const dvector& energies)
{
	double E1 = energies[imageNum - 1];
	double E2 = energies[imageNum];
	double E3 = energies[imageNum + 1];
	dvector tau(2 * L, 0);

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

	double length = std::inner_product(tau.begin(), tau.end(), tau.begin(), 0.0);
	for(int i = 0 ; i < 2 * L; ++i)
		tau[i] /= std::sqrt(length);
	return tau;

}

dvector calculateStepQuick(const dvector& perpV, int n, mat& velocities)
{
	dvector s(2 * L, 0);

	double fv = 0;
	double fd = 0;
	double mass = 1;
	double dt = 0.1;

	double tmp = 0;
	for(size_t i = 0; i < 2 * L; i++)
	{
		if( i == 0 || i == L )
			continue;    
		tmp = perpV[i] * velocities(n, i);
		if(tmp < 0) 
			velocities(n, i) = 0;
		else 
			fv += tmp;
		fd += perpV[i] * perpV[i];
	}

	tmp = 0;
	for(size_t i = 0; i < 2 * L; i++)
	{
		if( i == 0 || i == L )
			continue;
		velocities(n, i) = perpV[i] * (fv / fd + dt / mass);
		s[i] = velocities(n, i) * dt;
	}
	double length = std::inner_product(s.begin(), s.end(), s.begin(), 0.0);
	if( length > 1 )
		for(int i = 0 ; i < 2 * L; ++i)
			s[i] /= std::sqrt(length);
	return s;
}

dvector calculateStepQuick2(const std::vector< dvector >& path, const dvector& perpV, int n, mat& velocities)
{
	dvector s(2 * L, 0);

	double fv = 0;
	double fd = 0;
	double mass = 1.0;
	double dt = 0.1;

	double tmp = 0;
	for(size_t i = 0; i < 2 * L; i++)
	{
		//                if( i == 0 || i == L )
		//                    continue;
		tmp = perpV[i] * velocities(n, i);
		if(tmp < 0) velocities(n, i) = 0;
		else fv += tmp;
		fd += perpV[i] * perpV[i];
	}

	tmp = 0;
	for(size_t i = 0; i < 2 * L; i++)
	{
		//            if( i == 0 || i == L )
		//                continue;    
		velocities(n, i) = perpV[i] * (fv / fd + dt / mass);
		s[i] = velocities(n, i) * dt;
	}
	double length = std::inner_product(s.begin(), s.end(), s.begin(), 0.0);
	if( length > 1 )
		for(int i = 0 ; i < 2 * L; ++i)
			s[i] /= std::sqrt(length);
	return s;
}


void buldMEP(dvector thetaAngles, dvector phiAngles, dvector magneticMoments, dvector electronsNumber, 
		dvector E0, dvector U0, dmatrix hopingIntegrals, int numberOfImages, dvector initialState, dvector finalState, std::string fileName )
{
	dmatrix pathN(numberOfImages);
	dmatrix pathM(numberOfImages);    
	dmatrix gradients(numberOfImages);
	dvector energies(numberOfImages);
	mat velocities(numberOfImages, 2 * L, fill::zeros);
	mat path(numberOfImages, 2 * L, fill::zeros);

	for(int i = 0; i < 2 * L; ++i)
	{
		if( i == 0 || i == L )
			continue;

		double dh = (finalState[i] - initialState[i]) / (numberOfImages - 1);
		for(int j = 0; j < numberOfImages; ++j)
			path(j, i) = initialState[i] + dh * j;
	}  

	std::stringstream _name;
	_name <<  "initialState_" << fileName << ".txt";
	recordPath( path, _name.str().c_str() );        

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

	int maxEnergyPoint = 1;
	unsigned int iteration = 0;
	do
	{
		//        ++iteration;
		//        if( ++iteration > 500 )
		//            break;
		maxEnergyPoint = 1;

		for(int i = 1; i < numberOfImages - 1; ++i)
		{
			buildSelfconsistentSolution(i, path, pathM, pathN, E0, U0,  hopingIntegrals, gradients, energies); 
			if( energies[maxEnergyPoint] < energies[i])
				maxEnergyPoint = i;
		}

		dvector resultForce(2 * L, 0);
		for(int i = 1; i < numberOfImages - 1; ++i)
		{
			dvector tauVector = pathTangent(i, path, energies);

			dvector springForce(2 * L, 0);
			double gradientProjection    = 0;
			double springForceProjection = 0;

			for( int j = 0; j < 2 * L; ++j )
			{
				if( j == 0 || j == L )
					continue;      
				springForce[j]          =  0.5 * ( path( i + 1, j) + path(i - 1, j) - 2 * path(i, j) );
				gradientProjection      += gradients[i][j] * tauVector[j];
				springForceProjection   += springForce[j] * tauVector[j];
			}

			dvector perpV(2 * L, 0);

			for( int j = 0; j < 2 * L; ++j )
			{
				if( j == 0 || j == L )
					continue;
				if( i == maxEnergyPoint )
					perpV[j] = -gradients[i][j] + 2 * tauVector[j] * gradientProjection;
				else
					perpV[j] = -( gradients[i][j] - tauVector[j] * gradientProjection) + springForceProjection * tauVector[j];  

				resultForce[j] += perpV[j];                
			}

			dvector shift = calculateStepQuick(perpV, i, velocities);
			for(int j = 0; j < 2 * L; ++j)
				path(i, j) += shift[j];          
		}
		double totalForce = std::inner_product(resultForce.begin(), resultForce.end(), resultForce.begin(), 0.0);
		cout << "Total force: " << totalForce << std::endl;

		if( std::sqrt( std::abs(totalForce) )  < 1e-7 )
			break;
	} while( true );

	cout << "Final path:\n";
	for(int i = 0; i < numberOfImages; ++i )
		cout << path(i, 1)   << " " <<  path(i, 2) << " " << 0 << std::endl;

	cout << "Max energy point: " <<  maxEnergyPoint << " " << path(maxEnergyPoint, 1)  << " " << path(maxEnergyPoint, 2) << "\n";

	_name.clear();
	_name.str(std::string());
	_name <<  "finalState_" << fileName << ".txt";
	recordPath( path, _name.str().c_str() );        
	recordPoint( path, pathN, pathM, energies, "saddle.txt", maxEnergyPoint );
}

void recordPath( const mat& path, const char* fileName)
{
	std::ofstream logFile;
	logFile.open( fileName, std::ofstream::out | std::ofstream::trunc );
	for(int i = 0; i < path.n_rows; ++i )
		logFile << path(i, 2) << " " << path(i, 1) << " " << 0 << std::endl;
}

void recordPoint( const mat& path, const dmatrix& pathN, const dmatrix& pathM, const dvector& energies, const char* fileName, const int& pointIndex )
{
	std::ofstream logFile;
	logFile.open( fileName, std::ofstream::out | std::ofstream::app );
	logFile << path(pointIndex, 2) << " " << path(pointIndex, 1) << " " << 0 <<  " ";
	for(int i = 0; i < pathM[pointIndex].size(); ++i)
		logFile << pathM[pointIndex][i] << " ";
	logFile << energies[pointIndex] << std::endl;
}

std::vector< dvector > findMinima(const dvector& magneticMoments, const dvector& electronsNumber, 
		const dvector& E0, const dvector& U0, const dmatrix& hopingIntegrals)
{
	const int dotNumber = 2;
	std::vector< dvector > dots( 0, dvector() );

	//    boost::random::mt19937 rng; 
	//    boost::random::uniform_real_distribution<> randAngle(0, 2 * boost::math::constants::pi<double>()); 
	//    for(int i = 0; i < 2 * L; ++i)
	//        dots[0].push_back( randAngle(rng) );
	//MN angles
	//    dots.push_back( {0, -1.81989, 2.18787, 0, 0 ,0});   
	//    dots.push_back( {0,  4.47895, 2.19924, 0, 0, 0} );
	//    dots.push_back( {0,  1.81986, 4.09529, 0, 0, 0} );
	//    dots.push_back( {0,  4.46324, 8.47101, 0, 0, 0} );
	//CR angles
	dots.push_back( {0,  3.0 * 3.14, 1.0 * 3.14, 0, 0, 0} );

	//FE angles
	std::ofstream logFile;
	logFile.open( "minimas.txt", std::ofstream::out | std::ofstream::trunc );

	for(int i = 0; i < dots.size(); ++i)
	{
		double totalForce;
		dvector N;
		dvector M; 
		dvector Energy;
		do
		{
			totalForce = 0;
			N = electronsNumber;
			M = magneticMoments; 
			Energy = E0;
			unsigned int iterations = 0;
			dvector thetaAngles(L, 0);
			dvector phiAngles(L, 0);
			for(int j = 0; j < L; ++j)
			{
				thetaAngles[j] = dots[i][j];
				phiAngles[j]   = dots[i][j + L];
			}        

			bool isConsist;
			do
			{
				isConsist = true;
				for (int j = 0; j < L; ++j)
					isConsist &= SConsist(j, thetaAngles, phiAngles, M, N, E0, U0, hopingIntegrals, Energy);
				if( ++iterations == 50 )
					throw "#Infinite consist cluster loop.";    
			} while( !isConsist );           
			dvector gradient = GradE(thetaAngles, phiAngles, M, N, E0, U0, hopingIntegrals);
			for(int j = 0; j < 2 * L; ++j)
			{
				if( j == 0 || j == L )
					continue;
				dots[i][j] += -1.5 * gradient[j];
				totalForce += gradient[j] * gradient[j];
			}
			std::cout << "Total force: " << totalForce << std::endl;
		} while( std::sqrt(totalForce) > 1e-4 );
		logFile << dots[i][0] << " " << dots[i][1] << " " << dots[i][2] << " " << M[0] << " " << M[1] << " " << M[2] << " " << std::accumulate( Energy.begin(), Energy.end(), 0.0 ) << std::endl;

	}
	//    std::ofstream logFile;
	//    logFile.open( "MN_min.txt", std::ofstream::out | std::ofstream::trunc );
	//    for(int i = 0; i < dots[0].size(); ++i )
	//            logFile << dots[0][i] << " ";
	//    for(int i = 0; i < dots.size(); ++i)
	//        std::cout << i << ")" << dots[i][1] << " " << dots[i][2] << "\n";

	return dots;
}

std::vector< dvector > findMinimaVelocities(const dvector& magneticMoments, const dvector& electronsNumber, 
		const dvector& E0, const dvector& U0, const dmatrix& hopingIntegrals)
{
	std::vector< dvector > dots( 0, dvector() );

	//MN angles
	//    dots.push_back( {0, -1.81, 2.1, 0, 0 ,0});   
	//    dots.push_back( {0,  4.47, 2.1, 0, 0, 0} );
	//    dots.push_back( {0,  1.81, 4.0, 0, 0, 0} );
	//    dots.push_back( {0,  4.46, 8.4, 0, 0, 0} );
	//CR angles
	dots.push_back( {0,  0.0       ,       3.14, 0, 0 ,0});   
	dots.push_back( {0,  0.5 * 3.14, 1.5 * 3.14, 0, 0, 0} );
	dots.push_back( {0,  1.5 * 3.14, 0.5 * 3.14, 0, 0, 0} );
	//FE angles
	//    dots.push_back( {0, -0.0 * 3.14, 0.0 * 3.14, 0, 0 ,0});   
	//    dots.push_back( {0,  2.0 * 3.14, 0.0 * 3.14, 0, 0 ,0});   
	//    dots.push_back( {0,  2.0 * 3.14, 2.0 * 3.14, 0, 0 ,0});   

	mat velocities(dots.size(), 2 * L, fill::zeros);
	dmatrix gradients(dots.size());

	std::ofstream logFile;
	logFile.open( "minimas.txt", std::ofstream::out | std::ofstream::trunc );

	for(int i = 0; i < dots.size(); ++i)
	{
		dvector N;
		dvector M; 
		dvector Energy;

		do
		{
			N = electronsNumber;
			M = magneticMoments; 
			Energy = E0;
			unsigned int iterations = 0;
			dvector thetaAngles(L, 0);
			dvector phiAngles(L, 0);
			for(int j = 0; j < L; ++j)
			{
				thetaAngles[j] = dots[i][j];
				phiAngles[j]   = dots[i][j + L];
			}        

			bool isConsist;
			do
			{
				isConsist = true;
				for (int j = 0; j < L; ++j)
					isConsist &= SConsist(j, thetaAngles, phiAngles, M, N, E0, U0, hopingIntegrals, Energy);
				if( ++iterations == 50 )
					throw "#Infinite consist cluster loop.";    
			} while( !isConsist );    

			dvector perpV(2 * L, 0);
			dvector gradient = GradE(thetaAngles, phiAngles, M, N, E0, U0, hopingIntegrals);
			double length = std::inner_product(gradient.begin(), gradient.end(), gradient.begin(), 0.0);
			if( std::sqrt(length) < 1e-8)
				break;     

			for(int j = 0; j < 2 * L; ++j)
			{
				if( j == 0 || j == L )
					continue;
				perpV[j] = -gradient[j];
			}
			dvector shift = calculateStepQuick2(dots, perpV, i, velocities);
			for(int j = 0; j < 2 * L; ++j)
				dots[i][j] += shift[j];
			//            std::cout << "Total force: " << length << " " << dots[i][1] << ":" << dots[i][2] << std::endl;
		} while( true );
		logFile << dots[i][0] << " " << dots[i][1] << " " << dots[i][2] << " " << M[0] << " " << M[1] << " " << M[2] << " " << std::accumulate( Energy.begin(), Energy.end(), 0.0 ) << std::endl;

	}
	return dots;
}

std::vector< dvector > findMinimaVelocitiesThreaded(const dvector& magneticMoments, const dvector& electronsNumber, 
		const dvector& E0, const dvector& U0, const dmatrix& hopingIntegrals)
{
	std::vector< dvector > dots( 0, dvector() );

	//MN angles
	//    dots.push_back( {0, -1.81, 2.1, 0, 0 ,0});   
	//    dots.push_back( {0,  4.47, 2.1, 0, 0, 0} );
	//    dots.push_back( {0,  1.81, 4.0, 0, 0, 0} );
	//    dots.push_back( {0,  4.46, 8.4, 0, 0, 0} )6;
	//CR angles
	dots.push_back( {0,  0   * boost::math::constants::pi<double>(), 1.0 * boost::math::constants::pi<double>(), 0, 0 ,0});   
	dots.push_back( {0,  0.5 * boost::math::constants::pi<double>(), 1.5 * boost::math::constants::pi<double>(), 0, 0, 0} );
	dots.push_back( {0,  1.5 * boost::math::constants::pi<double>(), 0.5 * boost::math::constants::pi<double>(), 0, 0, 0} );
	//FE angles
	//    dots.push_back( {0, -0.0 * 3.14, 0.0 * 3.14, 0, 0 ,0});   
	//    dots.push_back( {0,  2.0 * 3.14, 0.0 * 3.14, 0, 0 ,0});   
	//    dots.push_back( {0,  2.0 * 3.14, 2.0 * 3.14, 0, 0 ,0});   

	mat velocities(dots.size(), 2 * L, fill::zeros);
	dmatrix gradients(dots.size());

	std::ofstream logFile;
	logFile.open( "minimas.txt", std::ofstream::out | std::ofstream::trunc );

	for(int i = 0; i < dots.size(); ++i)
	{
		dvector N;
		dvector M; 
		dvector Energy;

		BOOST_LOG_TRIVIAL(info) << "Optimizing dot #" << i;
		do
		{
			N = electronsNumber;
			M = magneticMoments; 
			Energy = E0;
			unsigned int iterations = 0;
			dvector thetaAngles(L, 0);
			dvector phiAngles(L, 0);
			for(int j = 0; j < L; ++j)
			{
				thetaAngles[j] = dots[i][j];
				phiAngles[j]   = dots[i][j + L];
			}        

			bool isConsist;
			/*
			   do
			   {
			   isConsist = true;
			   for (int j = 0; j < L; ++j)
			   isConsist &= SConsist(j, thetaAngles, phiAngles, M, N, E0, U0, hopingIntegrals, Energy);
			   if( ++iterations == 50 )
			   throw "#Infinite consist cluster loop.";    
			   } while( !isConsist );
			   */
			do
			{
				isConsist = true;
				int atom = 0;
				while( atom < L )
				{	  
					boost::thread_group threads;

					std::pair< dvector, dvector > resultNM = std::make_pair(N, M);
					for ( int threadID = 0; threadID < CORE_COUNT && atom < L; ++threadID, ++atom )
					{
						threads.add_thread( new boost::thread( &SConsistThreaded, atom, boost::ref(resultNM), std::make_pair(thetaAngles, phiAngles), 
									M, N, std::make_pair(E0, U0), 
									hopingIntegrals, boost::ref(Energy), boost::ref(isConsist)  ) );
					}
					threads.join_all();
					if( interrupt == true )
					{
						interrupt = false;
						throw "#Infinite selfconsist procedure.";  
					}
					N = resultNM.first;
					M = resultNM.second;
				}
				if( ++iterations == 500 )
					throw "#Infinite consist cluster loop.";    
			} while (!isConsist);


			dvector perpV(2 * L, 0);
			dvector gradient = GradE(thetaAngles, phiAngles, M, N, E0, U0, hopingIntegrals);
			double length = std::inner_product(gradient.begin(), gradient.end(), gradient.begin(), 0.0);
			if( std::sqrt(length) < 1e-8)
				break;     

			for(int j = 0; j < 2 * L; ++j)
			{
				if( j == 0 || j == L )
					continue;
				perpV[j] = -gradient[j];
			}
			dvector shift = calculateStepQuick2(dots, perpV, i, velocities);
			for(int j = 0; j < 2 * L; ++j)
				dots[i][j] += shift[j];
			//            std::cout << "Total force: " << length << " " << dots[i][1] << ":" << dots[i][2] << std::endl;
		} while( true );
		logFile << dots[i][0] << " " << dots[i][1] << " " << dots[i][2] << " " << M[0] << " " << M[1] << " " << M[2] << " " << std::accumulate( Energy.begin(), Energy.end(), 0.0 ) << std::endl;

	}
	return dots;
}

namespace JobFunctions 
{

	void printNeighbours( const dmatrix& arr )
	{
		for(int i = 0; i < arr.size(); ++i)
		{
			std::cout << i + 1 << ": ";
			for(int j = 0; j < arr[i].size(); ++j)
				if( arr[i][j] != 0 )
					std::cout << j + 1 << " ";
			std::cout << std::endl;
		}
	}
}

void init()
{
	logging::core::get()->set_filter
		(
		 logging::trivial::severity >= logging::trivial::info
		);
}



int main(int argc, char* argv[])
{
	init();

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
		remove("saddle.txt");
		std::ifstream hparams_file("hopingTrimer.txt");
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

		dvector Energy(L, 0);
		bool isConsist = true;
		thetaAngles[1] = -3.141592653589793;
		thetaAngles[2] = -2.8530790425458323;
		unsigned int iterations = 0;
		do
		{
			isConsist = true;
			for (int i = 0; i < L; ++i)
				isConsist &= SConsist(i, thetaAngles, phiAngles, MagneticMoments, ElectronsNumber, E0, U0, hoping_integrals, Energy);
			if( ++iterations == 50 )
				throw "#Infinite consist cluster loop.";    
		}
		while (!isConsist);
		return 0;

		//        dmatrix results = buildEnergySurface(thetaAngles, phiAngles, MagneticMoments, ElectronsNumber, E0, U0, hoping_integrals, true );
		// 
		//        std::ofstream logFile;
		//        logFile.open( "temp.txt", std::ofstream::out | std::ofstream::trunc );
		//        for( int i = 0; i < results.size(); ++i )
		//        {
		//            for( int j = 0; j < results[i].size(); ++j )
		//                logFile << std::fixed << std::setprecision(10) << results[i][j] << " ";
		//            logFile << std::endl;
		//        }
		//        return 0;

		std::cout << "Finding minimas\n";

		//        std::vector< dvector > dots = findMinimaVelocities(MagneticMoments, ElectronsNumber, E0, U0, hoping_integrals);   

		//        std::cout << "dot1 th3 " << hoping_integrals[0][1] << " " << dots[0][2] << "\n";
		//        std::cout << "dot1 th2 " << hoping_integrals[0][1] << " " << dots[0][1] << "\n";
		//        std::cout << "dot2 th3 " << hoping_integrals[0][1] << " " << dots[1][2] << "\n";
		//        std::cout << "dot2 th2 " << hoping_integrals[0][1] << " " << dots[1][1] << "\n";
		//        std::cout << "dot3 th3 " << hoping_integrals[0][1] << " " << dots[2][2] << "\n";
		//        std::cout << "dot3 th2 " << hoping_integrals[0][1] << " " << dots[2][1] << "\n";

		for(int i = 0; i <= 0; ++i)
		{
			double buf = i * 0.05 + 0.6;
			hoping_integrals[0][1] = hoping_integrals[0][2] = buf;
			std::vector< dvector > dots = findMinimaVelocitiesThreaded(MagneticMoments, ElectronsNumber, E0, U0, hoping_integrals);   
			std::cout << hoping_integrals[0][1] << " " << dots[1][1] << "\n";
		}


		//        std::cout << "Building MEP\n";
		//        buldMEP(thetaAngles, phiAngles, MagneticMoments, ElectronsNumber, E0, U0, hoping_integrals, 15,
		//                        dots[1], dots[2], "second" );
		//        buldMEP(thetaAngles, phiAngles, MagneticMoments, ElectronsNumber, E0, U0, hoping_integrals, 15,
		//                        dots[0], dots[2], "first" );
		//        buldMEP(thetaAngles, phiAngles, MagneticMoments, ElectronsNumber, E0, U0, hoping_integrals, 15,
		//                        dots[2], dots[3], "third" );


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
