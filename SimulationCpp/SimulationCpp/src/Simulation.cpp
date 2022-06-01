#include <cstdlib> //for random
#include <ctime> //for random
#include <initializer_list>
#include <iostream>

#include "Util.h"
#include "Simulation.h"
#include "Particle3DSGP.h"
#include "Particle3DPWG.h"
#include "Sphere.h"
#include "Ellipsoid.h"
#include "Planes.h"
#include "CylinderBasic.h"

Simulation::Simulation(int nStep, double timeStep, double envT2, double envDiffusivity, double envSizeX, double envSizeY, double envSizeZ)
	: m_nStep{ nStep },
	m_timeStep{ timeStep },
	m_timeTolerance{timeStep*Simulation::m_TOL}, //time used for particles to move right after collision
	m_environment{ Environment{envT2, envDiffusivity, envSizeX, envSizeY, envSizeZ} }
{
	std::srand(static_cast<unsigned int>(std::time(nullptr)));
	std::rand();
}

void Simulation::addSphere(double x, double y, double z, double T2, double diffusivity, double permeability, double radius) {
	m_compartments.push_back(std::make_shared<Sphere>(Sphere{ x, y, z, T2, diffusivity, permeability, radius }));
}

void Simulation::addEllipsoid(double x, double y, double z, double T2, double diffusivity, double permeability,
	double a[3], int len_a, double b[3], int len_b, double c[3], int len_c) {
	m_compartments.push_back(std::make_shared<Ellipsoid>(Ellipsoid{ x, y, z, T2, diffusivity, permeability,
		Eigen::Vector3d{a},  Eigen::Vector3d{b}, Eigen::Vector3d{c} }));
}

void Simulation::addEllipsoid(double x, double y, double z, double T2, double diffusivity, double permeability,
	double a, double b, double c) {
	m_compartments.push_back(std::make_shared<Ellipsoid>(Ellipsoid{ x, y, z, T2, diffusivity, permeability,
		Eigen::Vector3d{a, 0, 0},  Eigen::Vector3d{0, b, 0}, Eigen::Vector3d{0, 0, c} }));
}

void Simulation::addPlanes(double T2, double diffusivity, double spacing) {
	m_compartments.push_back(std::make_shared<Planes>(Planes{ T2, diffusivity, spacing }));
}

void Simulation::addCylinderBasic(double x, double y, double T2, double diffusivity, double permeability, double radius) {
	m_compartments.push_back(std::make_shared<CylinderBasic>(CylinderBasic{ x, y, T2, diffusivity, permeability, radius }));
}

void Simulation::createStartingPositions(int nPart, bool insideCompartments) {
	//must be called after all compartments were added
	m_nPart = nPart;
	Eigen::Vector3d envSize = m_environment.getSize();
	if (insideCompartments) {
		//only adding particles inside the compartments
		int i = 0;
		while (i < nPart) {
			Eigen::Vector3d pos = Util::getRandomUniformPosition(envSize);
			for (std::shared_ptr<AbstractCompartment> comp : m_compartments) {
				if (comp->contains(pos)) {
					m_startingPositions.push_back(pos);
					i++;
					break;
				}
			}
		}
	}
	else {
		for (int i = 0; i < nPart; i++) {
			m_startingPositions.push_back(Util::getRandomUniformPosition(envSize));
		}
	}
	m_compartments.push_back(std::make_shared<Environment>(m_environment)); //adding environment as final compartment
	m_nComp = (int) size(m_compartments);
}

void Simulation::createStartingPositionsAtPos(int nPart, double x, double y, double z) {
	//must be called after all compartments were added
	m_nPart = nPart;
	for (int i = 0; i < nPart; i++ ) {
		m_startingPositions.push_back(Eigen::Vector3d{ x, y, z });
	}
	m_compartments.push_back(std::make_shared<Environment>(m_environment)); //adding environment as final compartment
	m_nComp = (int) size(m_compartments);
}

void Simulation::createSequenceSGP() {
	for (int i = 0; i < m_nPart; i++) {
		Eigen::Vector3d pos = m_startingPositions[i];
		m_particles.push_back(std::make_shared<Particle3DSGP>(Particle3DSGP{ pos[0], pos[1], pos[2] }));
	}
}

void Simulation::createSequencePWG(double* sequenceTimes, int len_sequenceTimes,
	const double *sequenceVectors, int rows_sequenceVectors, int cols_sequenceVectors) {
	//std::cout << "DIM 1 " << rows_sequenceVectors << std::endl;
	//std::cout << "1st item " << sequenceVectors[4] << std::endl;
	std::vector<Eigen::Vector3d> sequenceValues;
	std::vector<double> seqTimes;
	for (int i = 0; i < rows_sequenceVectors; i++) {
		sequenceValues.push_back(Eigen::Vector3d{ sequenceVectors[3 * i], sequenceVectors[3 * i + 1], sequenceVectors[3 * i + 2] });
		seqTimes.push_back(sequenceTimes[i]);
		//std::cout << "sequence value " << sequenceValues[i] << std::endl;
	}
	//std::cout << "Times " << seqTimes[0] << seqTimes[1] << seqTimes[2] << seqTimes[3] << std::endl;
	for (int i = 0; i < m_nPart; i++) {
		Eigen::Vector3d pos = m_startingPositions[i];
		m_particles.push_back(std::make_shared<Particle3DPWG>(Particle3DPWG{ pos[0], pos[1], pos[2],
			std::make_shared<std::vector<double>>(seqTimes), std::make_shared<std::vector<Eigen::Vector3d>>(sequenceValues) }));
	}
}

void Simulation::getResults(double** output_vector1, int* output_length1) {
	int resSize = (int) size(m_results);
	*output_length1 = resSize;
	*output_vector1 = new double[resSize];
	for (int i = 0; i < resSize; i++) {
		(*output_vector1)[i] = m_results[i];
	}
	//*output_vector1 = &m_results[0];
}

std::shared_ptr<AbstractCompartment> Simulation::findCompartment(Eigen::Vector3d pos, std::shared_ptr<AbstractCompartment> excludedComp) {
	for (std::shared_ptr<AbstractCompartment> compartment : m_compartments) {
		if (compartment != excludedComp && compartment->contains(pos)) {
			return compartment;
		}
	}
	std::cout << "Error: no parent compartment found!" << std::endl;
	return nullptr; //if nothing is found (should not happen)
}

void Simulation::run(int seed, int partPrintNumber) {
	//setup
	if (seed >= 0) {
		std::srand(seed);
	}
	
	for (int i = 0; i < m_nComp - 1; i++) {
		std::shared_ptr<AbstractCompartment> compartment = m_compartments[i];
		std::shared_ptr<AbstractCompartment> parent = findCompartment(compartment->getPos(), compartment);
		compartment->setParentComp(parent);
	}

	//run
	for (int p = 0; p < m_nPart; p++) {
		std::shared_ptr<AbstractParticle3D> particle = m_particles[p];
		particle->setVelocity(Util::getRandomDirection());
		std::shared_ptr<AbstractCompartment> newComp = findCompartment(particle->getPos());
		particle->changeCompartment(newComp, m_timeStep);
		for (int n = 0; n < m_nStep; n++) {
			nextStep(particle);
			if (partPrintNumber > 0 && (p + 1) % partPrintNumber == 0) {
				std::cout << "Particle " << p << std::endl; //TODO ADD TIME TO PRINT
			}
		}

		//results
		std::vector<double> particleResults = particle->getResults();
		m_results.insert(m_results.end(), particleResults.begin(), particleResults.end());
	}
}

void Simulation::nextStep(std::shared_ptr<AbstractParticle3D> particle) {
	double t = 0; //time elapsed during the step
	particle->setVelocity(particle->getSpeed()*Util::getRandomDirection()); //random direction at each step
	while (t < m_timeStep) {
		std::array<Eigen::Vector3d, 2> ray = { particle->getPos(), particle->getVelocity() / particle->getSpeed() };
		std::vector<double> distances(m_nComp);
		double maxDistance = particle->getSpeed()*(m_timeStep - t);
		for (int c = 0; c < m_nComp; c++) {
			distances[c] = m_compartments[c]->findIntersection(ray, maxDistance);
		}
		int firstIndex = std::min_element(distances.begin(), distances.end()) - distances.begin();
		std::shared_ptr<AbstractCompartment> compartment = m_compartments[firstIndex];
		double reachTime = distances[firstIndex] / particle->getSpeed();

		if (t + reachTime + m_timeTolerance < m_timeStep) {
			//if there is an intersection
			Eigen::Vector3d oldPos = particle->getPos(); //making a copy
			particle->move(reachTime);
			if (compartment->collide(particle, oldPos, ray[0] + distances[firstIndex] * ray[1], m_timeStep)) {
				//if there has been a teleportation of the particle
				std::shared_ptr<AbstractCompartment> newComp = findCompartment(particle->getPos());
				particle->changeCompartment(newComp, m_timeStep);
			}
			particle->move(m_timeTolerance);
			t += reachTime + m_timeTolerance;
		} 
		else {
			//if there is no intersection
			particle->move(m_timeStep - t);
			t = m_timeStep;
		}
	}
}
