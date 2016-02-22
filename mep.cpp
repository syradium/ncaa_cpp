//
// Created by user on 6/28/15.
//

#include "mep.h"


extern bool interrupt;


std::vector<double> pathTangent(const int& imageNum, const arma::mat& path, const std::vector<double> & energies)
{
	double E1 = energies[imageNum - 1];
	double E2 = energies[imageNum];
	double E3 = energies[imageNum + 1];
	std::vector<double> tau(2 * L, 0);

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
		tau[i] /= sqrt(length);
	return tau;

}

std::vector<double> calculateStepQuick(const std::vector<double> & perpV, int n, arma::mat& velocities)
{
	std::vector<double> s(2 * L, 0);

	double fv = 0;
	double fd = 0;
	double mass = 1;
	double dt = 0.1;

	double tmp = 0;
	for(std::size_t i = 0; i < 2 * L; i++)
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
	for(std::size_t i = 0; i < 2 * L; i++)
	{
		if( i == 0 || i == L )
			continue;
		velocities(n, i) = perpV[i] * (fv / fd + dt / mass);
		s[i] = velocities(n, i) * dt;
	}
	double length = std::inner_product(s.begin(), s.end(), s.begin(), 0.0);
	if( length > 1 )
		for(int i = 0 ; i < 2 * L; ++i)
			s[i] /= sqrt(length);
	return s;
}

std::vector<double> calculateStepQuick2(const std::vector<std::vector<double> >& path, const std::vector<double> & perpV, int n, arma::mat& velocities)
{
	std::vector<double> s(2 * L, 0);

	double fv = 0;
	double fd = 0;
	double mass = 1.0;
	double dt = 0.1;

	double tmp = 0;
	for(std::size_t i = 0; i < 2 * L; i++)
	{
		//                if( i == 0 || i == L )
		//                    continue;
		tmp = perpV[i] * velocities(n, i);
		if(tmp < 0) velocities(n, i) = 0;
		else fv += tmp;
		fd += perpV[i] * perpV[i];
	}

	tmp = 0;
	for(std::size_t i = 0; i < 2 * L; i++)
	{
		//            if( i == 0 || i == L )
		//                continue;
		velocities(n, i) = perpV[i] * (fv / fd + dt / mass);
		s[i] = velocities(n, i) * dt;
	}
	double length = std::inner_product(s.begin(), s.end(), s.begin(), 0.0);
	if( length > 1 )
		for(int i = 0 ; i < 2 * L; ++i)
			s[i] /= sqrt(length);
	return s;
}

void buildMEP(
		std::vector<double> thetaAngles, std::vector<double> phiAngles, std::vector<double> magneticMoments, std::vector<double> electronsNumber,
              std::vector<double> E0, std::vector<double> U0, std::vector<std::vector<double>> hopingIntegrals, int numberOfImages, std::vector<double> initialState,
              std::vector<double> finalState, std::string fileName)
{
	std::vector<std::vector<double>> pathN(numberOfImages);
	std::vector<std::vector<double>> pathM(numberOfImages);
	std::vector<std::vector<double>> gradients(numberOfImages);
	std::vector<double> energies(numberOfImages);
	arma::mat velocities(numberOfImages, 2 * L, arma::fill::zeros);
	arma::mat path(numberOfImages, 2 * L, arma::fill::zeros);

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

    BuildSelfconsistentSolution(0, path, pathM, pathN, E0, U0, hopingIntegrals, gradients, energies);
    BuildSelfconsistentSolution(numberOfImages - 1, path, pathM, pathN, E0, U0, hopingIntegrals, gradients, energies);

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
            BuildSelfconsistentSolution(i, path, pathM, pathN, E0, U0, hopingIntegrals, gradients, energies);
			if( energies[maxEnergyPoint] < energies[i])
				maxEnergyPoint = i;
		}

		std::vector<double> resultForce(2 * L, 0);
		for(int i = 1; i < numberOfImages - 1; ++i)
		{
			std::vector<double> tauVector = pathTangent(i, path, energies);

			std::vector<double> springForce(2 * L, 0);
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

			std::vector<double> perpV(2 * L, 0);

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

			std::vector<double> shift = calculateStepQuick(perpV, i, velocities);
			for(int j = 0; j < 2 * L; ++j)
				path(i, j) += shift[j];
		}
		double totalForce = std::inner_product(resultForce.begin(), resultForce.end(), resultForce.begin(), 0.0);
		std::cout << "Total force: " << totalForce << std::endl;

		if( sqrt( std::abs(totalForce) )  < 1e-7 )
			break;
	} while( true );

	std::cout << "Final path:\n";
	for(int i = 0; i < numberOfImages; ++i )
		std::cout << path(i, 1)   << " " <<  path(i, 2) << " " << 0 << std::endl;

	std::cout << "Max energy point: " <<  maxEnergyPoint << " " << path(maxEnergyPoint, 1)  << " " << path(maxEnergyPoint, 2) << "\n";

	_name.clear();
	_name.str(std::string());
	_name <<  "finalState_" << fileName << ".txt";
	recordPath( path, _name.str().c_str() );
	recordPoint( path, pathN, pathM, energies, "saddle.txt", maxEnergyPoint );
}

void recordPath( const arma::mat& path, const char* fileName)
{
	std::ofstream logFile;
	logFile.open( fileName, std::ios_base::out | std::ios_base::trunc );
	for(int i = 0; i < path.n_rows; ++i )
		logFile << path(i, 2) << " " << path(i, 1) << " " << 0 << std::endl;
}

void recordPoint( const arma::mat& path, const std::vector<std::vector<double>> & pathN, const std::vector<std::vector<double>> & pathM, const std::vector<double> & energies, const char* fileName, const int& pointIndex )
{
	std::ofstream logFile;
	logFile.open( fileName, std::ios_base::out | std::ios_base::app );
	logFile << path(pointIndex, 2) << " " << path(pointIndex, 1) << " " << 0 <<  " ";
	for(int i = 0; i < pathM[pointIndex].size(); ++i)
		logFile << pathM[pointIndex][i] << " ";
	logFile << energies[pointIndex] << std::endl;
}

std::vector<std::vector<double> > findMinima(const std::vector<double> & magneticMoments, const std::vector<double> & electronsNumber,
		const std::vector<double> & E0, const std::vector<double> & U0, const std::vector<std::vector<double>> & hopingIntegrals)
{
	const int dotNumber = 2;
	std::vector<std::vector<double> > dots( 0, std::vector<double>() );

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
	logFile.open( "minimas.txt", std::ios_base::out | std::ios_base::trunc );

	for(int i = 0; i < dots.size(); ++i)
	{
		double totalForce;
		std::vector<double> N;
		std::vector<double> M;
		std::vector<double> Energy;
		do
		{
			totalForce = 0;
			N = electronsNumber;
			M = magneticMoments;
			Energy = E0;
			unsigned int iterations = 0;
			std::vector<double> thetaAngles(L, 0);
			std::vector<double> phiAngles(L, 0);
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
			std::vector<double> gradient = GradE(thetaAngles, phiAngles, M, N, E0, U0, hopingIntegrals);
			for(int j = 0; j < 2 * L; ++j)
			{
				if( j == 0 || j == L )
					continue;
				dots[i][j] += -1.5 * gradient[j];
				totalForce += gradient[j] * gradient[j];
			}
			std::cout << "Total force: " << totalForce << std::endl;
		} while( sqrt(totalForce) > 1e-4 );
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

std::vector<std::vector<double> > findMinimaVelocities(const std::vector<double> & magneticMoments, const std::vector<double> & electronsNumber,
		const std::vector<double> & E0, const std::vector<double> & U0, const std::vector<std::vector<double>> & hopingIntegrals)
{
	std::vector<std::vector<double> > dots( 0, std::vector<double>() );

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

	arma::mat velocities(dots.size(), 2 * L, arma::fill::zeros);
	std::vector<std::vector<double>> gradients(dots.size());

	std::ofstream logFile;
	logFile.open( "minimas.txt", std::ios_base::out | std::ios_base::trunc );

	for(int i = 0; i < dots.size(); ++i)
	{
		std::vector<double> N;
		std::vector<double> M;
		std::vector<double> Energy;

		do
		{
			N = electronsNumber;
			M = magneticMoments;
			Energy = E0;
			unsigned int iterations = 0;
			std::vector<double> thetaAngles(L, 0);
			std::vector<double> phiAngles(L, 0);
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

			std::vector<double> perpV(2 * L, 0);
			std::vector<double> gradient = GradE(thetaAngles, phiAngles, M, N, E0, U0, hopingIntegrals);
			double length = std::inner_product(gradient.begin(), gradient.end(), gradient.begin(), 0.0);
			if( sqrt(length) < 1e-8)
				break;

			for(int j = 0; j < 2 * L; ++j)
			{
				if( j == 0 || j == L )
					continue;
				perpV[j] = -gradient[j];
			}
			std::vector<double> shift = calculateStepQuick2(dots, perpV, i, velocities);
			for(int j = 0; j < 2 * L; ++j)
				dots[i][j] += shift[j];
			//            std::cout << "Total force: " << length << " " << dots[i][1] << ":" << dots[i][2] << std::endl;
		} while( true );
		logFile << dots[i][0] << " " << dots[i][1] << " " << dots[i][2] << " " << M[0] << " " << M[1] << " " << M[2] << " " << std::accumulate( Energy.begin(), Energy.end(), 0.0 ) << std::endl;

	}
	return dots;
}

std::vector<std::vector<double> > findMinimaVelocitiesThreaded(const std::vector<double> & magneticMoments, const std::vector<double> & electronsNumber,
		const std::vector<double> & E0, const std::vector<double> & U0, const std::vector<std::vector<double>> & hopingIntegrals)
{
	std::vector<std::vector<double> > dots( 0, std::vector<double>() );

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

	arma::mat velocities(dots.size(), 2 * L, arma::fill::zeros);
	std::vector<std::vector<double>> gradients(dots.size());

	std::ofstream logFile;
	logFile.open( "minimas.txt", std::ios_base::out | std::ios_base::trunc );

	for(int i = 0; i < dots.size(); ++i)
	{
		std::vector<double> N;
		std::vector<double> M;
		std::vector<double> Energy;

        cout << "Optimizing dot #" << i << std::endl;
		do
		{
			N = electronsNumber;
			M = magneticMoments;
			Energy = E0;
			unsigned int iterations = 0;
			std::vector<double> thetaAngles(L, 0);
			std::vector<double> phiAngles(L, 0);
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

					std::pair<std::vector<double>, std::vector<double> > resultNM = std::make_pair(N, M);
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


			std::vector<double> perpV(2 * L, 0);
			std::vector<double> gradient = GradE(thetaAngles, phiAngles, M, N, E0, U0, hopingIntegrals);
			double length = std::inner_product(gradient.begin(), gradient.end(), gradient.begin(), 0.0);
			if( sqrt(length) < 1e-8)
				break;

			for(int j = 0; j < 2 * L; ++j)
			{
				if( j == 0 || j == L )
					continue;
				perpV[j] = -gradient[j];
			}
			std::vector<double> shift = calculateStepQuick2(dots, perpV, i, velocities);
			for(int j = 0; j < 2 * L; ++j)
				dots[i][j] += shift[j];
			//            std::cout << "Total force: " << length << " " << dots[i][1] << ":" << dots[i][2] << std::endl;
		} while( true );
		logFile << dots[i][0] << " " << dots[i][1] << " " << dots[i][2] << " " << M[0] << " " << M[1] << " " << M[2] << " " << std::accumulate( Energy.begin(), Energy.end(), 0.0 ) << std::endl;

	}
	return dots;
}