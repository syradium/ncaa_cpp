//
// Created by user on 6/30/15.
//

#include "misc.h"

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