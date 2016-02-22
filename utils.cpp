#include "utils.h"


void GetParam(std::ifstream& file, std::vector<double>& arr)
{
    std::string buf;
    getline(file, buf);
    std::istringstream line(buf);
    copy(std::istream_iterator<double>(line), std::istream_iterator<double>(), back_inserter(arr));
}

void Init(dmatrix &hoping_integrals, dvector &E0, dvector &U0, dvector &ElectronsNumber, dvector &MagneticMoments,
		  dvector &thetaAngles, dvector &phiAngles) {
	remove("assets/saddle.txt");
	std::ifstream hparams_file("assets/hopingTrimer.txt");
	std::string buf;
	while (getline(hparams_file, buf))
	{
		std::istringstream line(buf);
		hoping_integrals.push_back(dvector());
		copy(std::istream_iterator<double>(line), std::istream_iterator<double>(), back_inserter(hoping_integrals.back()));
	}
	//Reading all necesary params from the text file
	std::ifstream init_file("assets/InitFile.txt");

	GetParam(init_file, E0);
	GetParam(init_file, U0);
	GetParam(init_file, ElectronsNumber);
	GetParam(init_file, MagneticMoments);
	GetParam(init_file, thetaAngles);
	GetParam(init_file, phiAngles);
}