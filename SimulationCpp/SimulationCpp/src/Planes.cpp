#include "Planes.h"
#include "Simulation.h"

Planes::Planes(double T2, double diffusivity, double spacing)
	: AbstractCompartment{0, 0, 0, T2, diffusivity},
	m_halfSpacing{spacing/2}
{
}

bool Planes::contains(Eigen::Vector3d pos) {
	return pos[0] > -m_halfSpacing && pos[0] < m_halfSpacing;
}

double Planes::findIntersection(std::array<Eigen::Vector3d, 2> ray, double maxDistance) {
	Eigen::Vector3d direction = ray[1];
	double dir_fraction = std::numeric_limits<double>::infinity();
	if (direction[0] != 0) {
		dir_fraction = 1 / direction[0];
	}
	
	double t1 = (m_halfSpacing - ray[0][0])*dir_fraction;
	double t2 = (-m_halfSpacing - ray[0][0])*dir_fraction;

	if (t1 > t2) { //MISTAKE FIXED???
		double temp = t1;
		t1 = t2;
		t2 = temp;
	}
	if (t1 > 0) {
		return t1;
	}
	else if (t2 > 0) {
		return t2;
	}
	return std::numeric_limits<double>::infinity();
}

bool Planes::collide(std::shared_ptr<AbstractParticle3D> particle, Eigen::Vector3d oldPos, Eigen::Vector3d intersection, double timeStep) {
	Eigen::Vector3d v = particle->getVelocity();
	v[0] = -v[0];
	particle->setVelocity(v);
	return false;
}
