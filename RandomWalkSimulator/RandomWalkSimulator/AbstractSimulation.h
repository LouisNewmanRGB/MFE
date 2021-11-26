#ifndef SIMULATION_H
#define SIMULATION_H

#include <Eigen/Dense>
#include <vector>

#include "AbstractParticle3D.h"
#include "Environment.h"
#include "AbstractCompartment.h"

class AbstractCompartment;

class AbstractSimulation
{
protected:
	std::vector<std::shared_ptr<AbstractParticle3D>> m_particles;
	int m_nPart;
	Environment m_environment;
	std::vector<std::shared_ptr<AbstractCompartment>> m_compartments;
	int m_nComp;
	int m_nStep;
	double m_timeStep;
	double m_timeTolerance;
	std::vector<Eigen::Vector3d> m_startingPositions;
	std::vector<double> m_results; //startingPosX, startingPosY, startingPosZ, posX, posY, posZ,
	//truePosX, truePosY, truePosz, t2signal, (sometimes) phase

public:
	static constexpr double m_TOL{ 1e-10 };

	AbstractSimulation(int nStep, double timeStep, std::vector<Eigen::Vector3d> startingPositions, Environment environment, std::vector<std::shared_ptr<AbstractCompartment>> compartments = {});

	std::shared_ptr<AbstractCompartment> findCompartment(Eigen::Vector3d pos, std::shared_ptr<AbstractCompartment> excludedComp = nullptr);
	void run(int seed = -1, int partPrintNumber = -1);
	void nextStep(std::shared_ptr<AbstractParticle3D> particle);

	virtual void generateParticles() = 0;
	virtual void saveResults(const std::string &path) = 0;

	std::vector<std::shared_ptr<AbstractParticle3D>> getParticles() { return m_particles; }
	std::vector<double> getResults() { return m_results; }
	double getTimeStep() { return m_timeStep; }

};

#endif