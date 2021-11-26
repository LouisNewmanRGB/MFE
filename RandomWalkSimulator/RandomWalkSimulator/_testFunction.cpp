#include "Scripts.h"

#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <cstdlib> //for random
#include <ctime>
#include <initializer_list>

#include "Environment.h"
#include "SimulationSGP.h"
#include "Sphere.h"
#include "Util.h"

void testFunction() {
	int nPart = 100000;
	int nStep = 8;
	double diffusionTime = 50; //ms
	double timeStep = diffusionTime / nStep;
	double envSize = 50; //um
	double D = 2.;
	double T2 = std::numeric_limits<double>::infinity();

	double radius = 8; //um
	double permeability = 0;

	std::vector<std::shared_ptr<AbstractCompartment>> compartments = { std::make_shared<Sphere>(Sphere{0, 0, 0, T2, D, permeability, radius}) };

	//std::shared_ptr<AbstractCompartment> env = std::make_shared<Environment>( T2, D, envSize, envSize, envSize );
	Environment env = Environment{ T2, D, envSize, envSize, envSize };
	//std::shared_ptr<AbstractCompartment> env2 = std::make_shared<Environment>(env);

	std::vector<Eigen::Vector3d> startingPos;

	std::srand(static_cast<unsigned int>(std::time(nullptr)));
	std::rand();

	//for (int k = 0; k < 10; k++) {
	//	std::cout << Util::getRandomUniform(-1, 2) << std::endl;
	//}

	for (int i = 0; i < nPart; i++) {
		//part.push_back( Particle3D{ Util::getRandomUniform(-envSize/2, envSize/2), Util::getRandomUniform(-envSize / 2, envSize / 2), Util::getRandomUniform(-envSize / 2, envSize / 2) } );
		startingPos.push_back(Util::getRandomQuadratic(radius)*Util::getRandomDirection());
		//part.push_back(Particle3D{ 0., 0., 0. });
	}

	SimulationSGP sim = SimulationSGP{ nStep, timeStep, startingPos, env , compartments };
	std::cout << "Done building simulation!" << std::endl;
	sim.run();
	std::cout << "Done simulating!" << std::endl;
	sim.saveResults("../results/npyTestSphere.npy");
	std::cout << "Done exporting!" << std::endl;
}