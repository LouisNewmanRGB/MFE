#include "Scripts.h"

#include <ctime>
#include <iostream>

#include "SimulationSGP.h"
#include "Planes.h"
#include "Util.h"

void planesMaxRadius() {
	int nPart = 1000000;
	int nStep = 8;
	double diffusionTime = 400; //ms
	double timeStep = diffusionTime / nStep;
	double D = 2.;
	double T2 = std::numeric_limits<double>::infinity();

	double spacing = 40; //um
	//double radius2 = radius * radius;
	double envSize = 5 * spacing; //um
	double permeability = 0;

	////////////////RANDOM SEEDING//////////////////
	std::srand(static_cast<unsigned int>(std::time(nullptr)));
	std::rand();
	////////////////////////////////////////////////

	std::vector<std::shared_ptr<AbstractCompartment>> compartments = { std::make_shared<Planes>(Planes{T2, D, spacing}) };
	Environment env = Environment{ T2, D, envSize, envSize, envSize };

	std::vector<Eigen::Vector3d> startingPos;
	for (int i = 0; i < nPart; i++) {
		startingPos.push_back(Eigen::Vector3d{Util::getRandomUniform(-spacing/2, spacing/2), Util::getRandomUniform(-envSize / 2, envSize / 2), Util::getRandomUniform(-envSize / 2, envSize / 2) } );
	}

	SimulationSGP sim = SimulationSGP{ nStep, timeStep, startingPos, env , compartments };
	std::cout << "Done building simulation!" << std::endl;
	sim.run();
	std::cout << "Done simulating!" << std::endl;
	sim.saveResults("../results/planes_max_radius.npy");
	std::cout << "Done exporting!" << std::endl;

}