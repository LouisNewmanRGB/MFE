#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "SimulationFiniteGradient.h"
#include "Particle3DPWG.h"
#include "npy.hpp"

SimulationFiniteGradient::SimulationFiniteGradient(int nStep, double timeStep, std::vector<Eigen::Vector3d> startingPositions, Environment environment,
	std::vector<std::shared_ptr<AbstractCompartment>> compartments, std::vector<double> sequenceTimes, std::vector<Eigen::Vector3d> sequenceValues)
	:AbstractSimulation{nStep, timeStep, startingPositions, environment, compartments},
	m_sequenceTimes{ sequenceTimes },
	m_sequenceValues{ sequenceValues }
{
}

void SimulationFiniteGradient::generateParticles() {
	for (int i = 0; i < m_nPart; i++) {
		Eigen::Vector3d pos = m_startingPositions[i];
		m_particles.push_back(std::make_shared<Particle3DPWG>(Particle3DPWG{ pos[0], pos[1], pos[2],
			std::make_shared<std::vector<double>>(m_sequenceTimes), std::make_shared<std::vector<Eigen::Vector3d>>(m_sequenceValues) }));
	}
}

void SimulationFiniteGradient::saveResults(const std::string &path) {
	const long unsigned shapeArray[] = { m_nPart, 11 };
	npy::SaveArrayAsNumpy(path, false, 2, shapeArray, m_results);
}