#include <cstdlib> //for random
#include <initializer_list>
#include <iostream>

#include "Util.h"
#include "AbstractSimulation.h"
#include "Particle3DSGP.h"

AbstractSimulation::AbstractSimulation(int nStep, double timeStep, std::vector<Eigen::Vector3d> startingPositions, Environment environment, std::vector<std::shared_ptr<AbstractCompartment>> compartments)
	: m_nStep{ nStep },
	m_timeStep{ timeStep },
	m_startingPositions{ startingPositions },
	m_nPart{ (int) size(startingPositions) },
	m_environment{ environment },
	m_nComp{ (int) size(compartments) + 1},
	m_timeTolerance{timeStep*AbstractSimulation::m_TOL} //time used for particles to move right after collision
{
	m_compartments = compartments;
	m_compartments.resize(m_nComp);
	m_compartments[m_nComp - 1] = std::make_shared<Environment>(environment);
}

std::shared_ptr<AbstractCompartment> AbstractSimulation::findCompartment(Eigen::Vector3d pos, std::shared_ptr<AbstractCompartment> excludedComp) {
	for (std::shared_ptr<AbstractCompartment> compartment : m_compartments) {
		if (compartment != excludedComp && compartment->contains(pos)) {
			return compartment;
		}
	}
	std::cout << "Error: no parent compartment found!" << std::endl;
	return nullptr; //if nothing is found (should not happen)
}

void AbstractSimulation::run(int seed, int partPrintNumber) {
	//setup
	if (seed >= 0) {
		std::srand(seed);
	}
	generateParticles();
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

void AbstractSimulation::nextStep(std::shared_ptr<AbstractParticle3D> particle) {
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
			compartment->collide(particle, oldPos, ray[0] + distances[firstIndex] * ray[1], *this);
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
