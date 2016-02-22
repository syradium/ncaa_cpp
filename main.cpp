/* 
 * File:   main.cpp
 * Author: Radium
 *
 * Created on January 14, 2013, 1:15 AM
 */

#include "definitions.h"
#include "gradient.h"
#include "mep.h"
#include "constants.h"
#include "threaded.h"

namespace logging = boost::log;

std::pair< dvector, dvector > ThreeDiag(const hermitian_matrix<complex<double>, upper>& H, int pIndex, int qIndex, const complex<double>& k, const complex<double>& m)
{
	vector< std::complex<double> > y0(2 * L, 0);
	vector< std::complex<double> > y1(2 * L, 0);
	vector< std::complex<double> > y(2 * L, 0);
	//Subdiagonal coefficients of matrix
	vector< std::complex<double> > a(2 * L, 0);
	//Diagonal coefficiencrs of matrix
	vector< std::complex<double> > b(2 * L - 1, 0);
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
//	vector< std::complex<double> > buf = a[0] * y0;
	//Caution: y is not actually a basis vector. Here it is a y1 * b
	y = Hy - (a[0] * y0);
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
		aa.push_back(i.real());

	BOOST_FOREACH(std::complex<double> i, b)
		bb.push_back(i.real() * i.real() );
	aa.resize(imax + 1);
	return std::make_pair(aa, bb);
}

double find_root(double a, double b, const double &a0, const double &b0, const vector<double> &p0,
				 const vector<double> &q0, const int &ii0)
{
	auto f = [&a0, &b0, &p0, &q0, &ii0](double x) -> double
	{
		double s = 0;
		for (int i = 0; i <= ii0; ++i)
			s += p0[i] / (x - q0[i]);
		return x - a0 - b0 * s;
	};

	auto result = bisect(f, a, b, [](double x, double y) -> bool { return y - x <= ::eps; });
	return result.first;
}

int FR(const double& a, const double& b, matrix<double>& t, vector<double>& p0, vector<double>& q0, const int& ii0)
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
		return 0;
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
				p0[j] += t(0, i );
		}
	}

	double leftBound = q0[0] - 1e4;
	double rightBound = q0[j] + 1e4;

	z[0] = find_root(leftBound, q0[0] - ::eps * 10, a0, b0, p0, q0, j);

	for (int i = 0; i < j; ++i)
		z[i + 1] = find_root(q0[i] + ::eps * 10, q0[i + 1] - ::eps * 10, a0, b0, p0, q0, j);
	z[j + 1] = find_root(q0[j] + ::eps * 10, rightBound, a0, b0, p0, q0, j);

	for (int i = 0; i <= j + 1; ++i)
	{
		double pTerm = 1;
		double pDenom = 1;
		for (int p = 0; p <= j; ++p)
		{
			pTerm *= z[i] - q0[p];
			if (i != p)
				pDenom *= z[i] - z[p];
		}
		if (i != j + 1)
			pDenom *= z[i] - z[j + 1];
		t(0, i) = pTerm / pDenom;
		t(1, i) = z[i];
	}
	return j + 1;
}


//Calculates the coefficients P and Q of the Green`s function`s matrix element represented by a fraction
//"a" and "b" are the diagonal and subdiagonal elemnts of the hamiltonian
matrix<double> Green(const dvector& a, const dvector& b)
{
	matrix<double> t(2, 2 * L, 0);
	vector<double> p0(2 * L, 0);
	vector<double> q0(2 * L, 0);

	t(0, 0) = 1;
	t(1, 0) = a.back();

	int ii00 = 0;
	for (int i = a.size() - 2; i >= 0; --i)
		ii00 = FR(a[i], b[i], t, p0, q0, ii00);
	t.resize(2, ii00 + 1);
	return t;
}

//Carry out an integration of the density matrix, which is the number of d-electrons, if a diagonal element is considered
double CN(const matrix<double>& t)
{
	double s1 = 0;
	double s2 = 0;
	for (int i = 0; i < t.size2(); ++i)
	{
		s1 += t(0, i);
		s2 += t(0, i) * std::atan(t(1, i));
	}

	return 0.5 - s2 * 0.318309886183790671;
}

double CalculateEnergy(const matrix<double>& t)
{
	double s1 = 0;
	double s2 = 0;
	double s3 = 0;

	for (int i = 0; i < t.size2(); ++i)
	{
		double PQ = t(0, i) * t(1, i);
		s1 += PQ;
		s2 += PQ * atan( t(1, i) );
		s3 += t(0, i) * log(t(1, i) * t(1, i) + 1);
	}
	return 0.5 * s1 - (s2 - 0.5 * s3) * 0.318309886183790671;
}

bool SConsist(int i, const dvector& tAngle, const dvector& pAngle, dvector& M, dvector& N, const dvector& E0, const dvector& U0, const dmatrix& hopingIntegrals, dvector& E)
{
	double magneticField = 0;
	//Saving magnetic moments for each atom
	double m_new = M[i];
	//Saving number of d-electrons on each atom
	double n_new = N[i];
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
		N[i] = n_new;
		M[i] = m_new;

		H.clear(); greenFraction.clear();

		formHamiltonian(H, tAngle, pAngle, N, M, E0, U0, hopingIntegrals, magneticField); 

		auto green_matrix_element = [&H](int index1, int index2, complex<double> vector1, complex<double> vector2) -> matrix<double>
		{
			auto cofficients = ThreeDiag(H, index1, index2, vector1, vector2);
			return Green(cofficients.first, cofficients.second);
		};

		greenFraction = green_matrix_element(i, i, complex<double>(1), complex<double>(1));
		double Nu = CN(greenFraction);
		double Eu = CalculateEnergy(greenFraction);

		greenFraction = green_matrix_element(i + L, i + L, complex<double>(1), complex<double>(1));
		double Nd = CN(greenFraction);
		double Ed = CalculateEnergy(greenFraction);

		greenFraction = green_matrix_element(i, i + L, complex<double>(1), complex<double>(1));
		double SFp = CN(greenFraction);

		greenFraction = green_matrix_element(i, i + L, complex<double>(1), complex<double>(-1));
		double SFn = CN(greenFraction);

		greenFraction = green_matrix_element(i, i + L, complex<double>(0, 1), complex<double>(1));
		double AFp = CN(greenFraction);

		greenFraction = green_matrix_element(i, i + L, complex<double>(0, 1), complex<double>(-1));
		double AFn = CN(greenFraction);

		n_new = Nu + Nd;
		m_new = (Nu - Nd) * cos(tAngle[i]) - ((SFp - SFn) * cos(pAngle[i]) - (AFp - AFn) * sin(pAngle[i])) * sin(tAngle[i]);
		E[i] = Eu + Ed - U0[i] * ( n_new * n_new - m_new * m_new ) * 0.25 ;
	}
	while ( (std::abs(m_new - M[i]) > delta || std::abs(n_new - N[i]) > delta)  && itr < maxIterCount ) ;
	N[i] = n_new;
	M[i] = m_new;

	if( itr == maxIterCount )
		throw "#Infinite selfconsist procedure.";
	return itr == 1;
}

dmatrix buildEnergySurface(dvector thetaAngles, dvector phiAngles, const dvector& magneticMoments, const dvector& electronsNumber,
		const dvector& E0, const dvector& U0, const dmatrix& hopingIntegrals, const bool& threaded )
{
	std::vector<dvector> results;
	int stepNumber          =  50;
	double theta2Begin      = -1.0 * boost::math::constants::pi<double>();
	double theta2End        =  2.0 * boost::math::constants::pi<double>();
	double theta3Begin      = -1.0 * boost::math::constants::pi<double>();
	double theta3End        =  2.0 * boost::math::constants::pi<double>();

	BOOST_LOG_TRIVIAL(info) << "Building energy surface with " << stepNumber << " steps";
	for(double th2: linspace(theta2Begin, theta2End, stepNumber))
	{
		BOOST_LOG_TRIVIAL(info) << "Step th2 " << th2 + 1 << " of " << stepNumber;

		dvector bufResults;
		thetaAngles[0] = 0;
		thetaAngles[1] = th2;

		for(double th3: linspace(theta3Begin, theta3End, stepNumber))
		{
			try
			{
				//Energy
				dvector Energy(L, 0);
				//initial electon number
				dvector N(electronsNumber);
				//initial magnenic moments
				dvector M(magneticMoments);

				thetaAngles[2] = th3;

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
//							if( interrupt == true )
//							{
//								interrupt = false;
//								throw "#Infinite selfconsist procedure.";
//							}
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

            }
			catch (const char* msg)
			{
				BOOST_LOG_TRIVIAL(error) << msg;
				std::terminate();
			}
			catch( boost::thread_interrupted const& e )
			{
                BOOST_LOG_TRIVIAL(error) << "Thread was interrupted";
				std::terminate();
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

	return results;
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

void GetParam(std::ifstream& file, std::vector<double>& arr)
{
    std::string buf;
    std::getline(file, buf);
    std::istringstream line(buf);
    std::copy(std::istream_iterator<double>(line), std::istream_iterator<double>(), std::back_inserter(arr));
}

void init()
{
	logging::core::get()->set_filter(logging::trivial::severity >= logging::trivial::info);
}


int main(int argc, char* argv[])
{
    std::cout << "Begin" << std::endl;
	init();
	std::cout << std::setprecision(16);
	//Hoping itegrals between atoms
	dmatrix hoping_integrals;
	//Position of non-perturbed d-level with respect to the Fermi energy for all atoms
	dvector E0;
	//Coulomb repulsion energy for all atoms
	dvector U0;
	//initial number of d-electrons
	dvector ElectronsNumber;
	//initial magnenic moments
	dvector MagneticMoments;
	//Theta angles for atoms
	dvector thetaAngles;
	//Phi angle for atoms
	dvector phiAngles;
	//Magnetic field
	double dh = 0;

	try 
	{
		remove("assets/saddle.txt");
		std::ifstream hparams_file("assets/hopingTrimer.txt");
		std::string buf;
		while (std::getline(hparams_file, buf))
		{
			std::istringstream line(buf);
			hoping_integrals.push_back(std::vector<double>());
			std::copy(std::istream_iterator<double>(line), std::istream_iterator<double>(), std::back_inserter(hoping_integrals.back()));
		}
		//Reading all necesary params from the text file
		std::ifstream init_file("assets/InitFile.txt");

		GetParam(init_file, E0);
		GetParam(init_file, U0);
		GetParam(init_file, ElectronsNumber);
		GetParam(init_file, MagneticMoments);
		GetParam(init_file, thetaAngles);
		GetParam(init_file, phiAngles);

		/*
		dvector Energy(L, 0);
		bool isConsist = true;
		thetaAngles[1] = -3.141592653589793;
		thetaAngles[2] = -3.141592653589793;
		unsigned int iterations = 0;
		do
		{
			isConsist = true;
			for (int i = 0; i < L; ++i)
			{
				bool result = SConsist(i, thetaAngles, phiAngles, MagneticMoments, ElectronsNumber, E0, U0, hoping_integrals, Energy);
				isConsist &= result;
			}
			if( ++iterations == 50 )
				throw "#Infinite consist cluster loop.";
		}
		while (!isConsist);
		return 0;
		*/

		        dmatrix results = buildEnergySurface(thetaAngles, phiAngles, MagneticMoments, ElectronsNumber, E0, U0, hoping_integrals, false );
		        std::ofstream logFile;
		        logFile.open( "temp.txt", std::ofstream::out | std::ofstream::trunc );
		        for( int i = 0; i < results.size(); ++i )
		        {
		            for( int j = 0; j < results[i].size(); ++j )
		                logFile << std::fixed << std::setprecision(10) << results[i][j] << " ";
		            logFile << std::endl;
		        }
		        return 0;

		std::cout << "Finding minimas\n";

		        std::vector< dvector > dots = findMinimaVelocities(MagneticMoments, ElectronsNumber, E0, U0, hoping_integrals);

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
		        buildMEP(thetaAngles, phiAngles, MagneticMoments, ElectronsNumber, E0, U0, hoping_integrals, 15,
		                        dots[1], dots[2], "second" );
		//        buildMEP(thetaAngles, phiAngles, MagneticMoments, ElectronsNumber, E0, U0, hoping_integrals, 15,
		//                        dots[0], dots[2], "first" );
		//        buildMEP(thetaAngles, phiAngles, MagneticMoments, ElectronsNumber, E0, U0, hoping_integrals, 15,
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
