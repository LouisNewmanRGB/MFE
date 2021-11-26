#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <npy.hpp> //for numpy saving

#include "SimulationSGP.h"
#include "Particle3DSGP.h"

SimulationSGP::SimulationSGP(int nStep, double timeStep, std::vector<Eigen::Vector3d> startingPositions, Environment environment,
	std::vector<std::shared_ptr<AbstractCompartment>> compartments)
	:AbstractSimulation(nStep, timeStep, startingPositions, environment, compartments)
{
}

void SimulationSGP::generateParticles() {
	//m_results.resize(m_nPart * 10); //TODO: try vec.reserve
	//m_particles.resize(m_nPart);
	for (int i = 0; i < m_nPart; i++) {
		Eigen::Vector3d pos = m_startingPositions[i];
		m_particles.push_back(std::make_shared<Particle3DSGP>(Particle3DSGP{ pos[0], pos[1], pos[2] }));
		//m_particles[i] = Particle3D{ pos[0], pos[1], pos[2] };
	}
}

void SimulationSGP::saveResults(const std::string &path) {
	const long unsigned shapeArray[] = { m_nPart, 10 };
	npy::SaveArrayAsNumpy(path, false, 2, shapeArray, m_results);
}