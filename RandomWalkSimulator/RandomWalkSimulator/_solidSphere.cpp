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

void solidSphere() {
	int nPart = 1000000;
	int nStep = 8;
	double diffusionTime = 5; //ms
	double timeStep = diffusionTime / nStep;
	double D = 2.;
	double T2 = std::numeric_limits<double>::infinity();

	double radius = 8; //um
	double radius2 = radius * radius;
	double envSize = 5*radius; //um
	double permeability = 0;

	////////////////RANDOM SEEDING//////////////////
	std::srand(static_cast<unsigned int>(std::time(nullptr)));
	std::rand();
	////////////////////////////////////////////////

	std::vector<std::shared_ptr<AbstractCompartment>> compartments = { std::make_shared<Sphere>(Sphere{0, 0, 0, T2, D, permeability, radius}) };
	Environment env = Environment{ T2, D, envSize, envSize, envSize };

	std::vector<Eigen::Vector3d> startingPos;

	/*
	for (int i = 0; i < nPart; i++) {
		//part.push_back( Particle3D{ Util::getRandomUniform(-envSize/2, envSize/2), Util::getRandomUniform(-envSize / 2, envSize / 2), Util::getRandomUniform(-envSize / 2, envSize / 2) } );
		//part.push_back(Particle3D{ 0., 0., 0. });
		startingPos.push_back(Util::getRandomQuadratic(radius)*Util::getRandomDirection());
	}
	*/
	int i = 0;
	while (i < nPart) {
		Eigen::Vector3d pos{ Util::getRandomUniform(-envSize / 2, envSize / 2), Util::getRandomUniform(-envSize / 2, envSize / 2), Util::getRandomUniform(-envSize / 2, envSize / 2) };
		if (pos.squaredNorm() > radius2) {
			startingPos.push_back(pos);
			i++;
		}
	}

	SimulationSGP sim = SimulationSGP{ nStep, timeStep, startingPos, env , compartments };
	std::cout << "Done building simulation!" << std::endl;
	sim.run();
	std::cout << "Done simulating!" << std::endl;
	sim.saveResults("../results/solid_sphere.npy");
	std::cout << "Done exporting!" << std::endl;
}