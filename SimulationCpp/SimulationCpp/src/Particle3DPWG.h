#ifndef PARTICLE3DPWG_H
#define PARTICLE3DPWG_H

#include "AbstractParticle3D.h"

#include <vector>

class Particle3DPWG : public AbstractParticle3D
{
private:
	double m_phase;
	double m_totalElapsedTime;
	int m_sequenceIndex;
	int m_maxIndex;
	double m_sequenceTime;
	Eigen::Vector3d m_sequenceValue;
	std::shared_ptr<std::vector<double>> m_sequenceTimes;
	std::shared_ptr<std::vector<Eigen::Vector3d>> m_sequenceValues;

	double integrateGradient(double time);

public:
	Particle3DPWG(double x, double y, double z, std::shared_ptr<std::vector<double>> sequenceTimes, std::shared_ptr<std::vector<Eigen::Vector3d>> sequenceValues);

	void move(double time);
	std::vector<double> getResults();
	
	double getPhase() { return m_phase; }

};

#endif

