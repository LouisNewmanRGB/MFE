#ifndef SIMULATION_H
#define SIMULATION_H

#include "Export.h"

#include <Eigen/Dense>
#include <vector>

#include "AbstractParticle3D.h"
#include "Environment.h"
#include "AbstractCompartment.h"

class AbstractCompartment;

class DLLEXPORT Simulation
{
protected:
	std::vector<std::shared_ptr<AbstractParticle3D>> m_particles;
	std::vector<Eigen::Vector3d> m_startingPositions;
	int m_nPart;
	Environment m_environment;
	std::vector<std::shared_ptr<AbstractCompartment>> m_compartments;
	int m_nComp;
	int m_nStep;
	double m_timeStep;
	double m_timeTolerance;
	std::vector<double> m_results; //startingPosX, startingPosY, startingPosZ, posX, posY, posZ,
	//truePosX, truePosY, truePosz, t2signal, (sometimes) phase

	void nextStep(std::shared_ptr<AbstractParticle3D> particle);
	std::shared_ptr<AbstractCompartment> findCompartment(Eigen::Vector3d pos, std::shared_ptr<AbstractCompartment> excludedComp = nullptr);

public:
	static constexpr double m_TOL{ 1e-10 };

	Simulation(int nStep, double timeStep, double envT2, double envDiffusivity, double envSizeX, double envSizeY, double envSizeZ);
	void addSphere(double x, double y, double z, double T2, double diffusivity, double permeability, double radius);
	void addEllipsoid(double x, double y, double z, double T2, double diffusivity, double permeability,
		double a[3], int len_a, double b[3], int len_b, double c[3], int len_c);
	void addEllipsoid(double x, double y, double z, double T2, double diffusivity, double permeability,
		double a, double b, double c);
	void addPlanes(double T2, double diffusivity, double spacing);
	void addCylinderBasic(double x, double y, double T2, double diffusivity, double permeability, double radius);
	void createStartingPositions(int nPart, bool insideCompartments);
	void createStartingPositionsAtPos(int nPart, double x, double y, double z);
	void createSequenceSGP();
	void createSequencePWG(double* sequenceTimes, int len_sequenceTimes, 
		const double *sequenceVectors, int rows_sequenceVectors, int cols_sequenceVectors);

	void run(int seed = -1, int partPrintNumber = -1);
	void getResults(double** output_vector1, int* output_length1); //returns output vector to Python

};

#endif