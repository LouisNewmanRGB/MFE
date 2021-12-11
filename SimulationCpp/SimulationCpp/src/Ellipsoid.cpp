#include "Ellipsoid.h"
#include "AbstractParticle3D.h"
#include "Util.h"
#include "Simulation.h"
#include <iostream>

Ellipsoid::Ellipsoid(double x, double y, double z, double T2, double diffusivity, double permeability, Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c)
	:AbstractCompartment{ x, y, z, T2, diffusivity },
	m_permeability{ permeability },
	m_maxRadius{ std::max({ a.norm(), b.norm(), c.norm() }) }
{
	Eigen::Matrix3d R;
	R.col(0) = a;
	R.col(1) = b;
	R.col(2) = c;
	m_invR = R.inverse();
	m_A = m_invR.transpose() * m_invR;
}

bool Ellipsoid::contains(Eigen::Vector3d pos) {
	return (m_invR*(pos - m_pos)).norm() < 1;
}

double Ellipsoid::findIntersection(std::array<Eigen::Vector3d, 2> ray, double maxDistance) {
	if ((ray[0] - m_pos).norm() - m_maxRadius <= maxDistance) {
		Eigen::Vector3d L = m_invR * (m_pos - ray[0]);
		Eigen::Vector3d invRD = m_invR * ray[1];
		double a = invRD.squaredNorm();
		double b = -2 * invRD.dot(L);
		double c = L.squaredNorm() - 1;
		double delta = b * b - 4 * a*c;
		if (delta > 0) {
			double t = (-b - sqrt(delta)) / (2 * a);
			if (t > 0) {
				return t;
			}
			else {
				double t = (-b + sqrt(delta)) / (2 * a);
				if (t > 0) {
					return t;
				}
			}
		}
		else if (delta == 0) {
			double t = -b / (2 * a);
			if (t > 0) {
				return t;
			}
		}
	}
	//no intersection
	return std::numeric_limits<double>::infinity();
}

bool Ellipsoid::collide(std::shared_ptr<AbstractParticle3D> particle, Eigen::Vector3d oldPos, Eigen::Vector3d intersection, double timeStep) {
	double transmissionProba = 4 * m_permeability / particle->getSpeed();
	if (transmissionProba > 1) {
		std::cout << "Error: greater than one transmission probability!" << std::endl;
	}
	assert(transmissionProba <= 1);

	if (Util::getRandomUniform(0, 1) >= transmissionProba) {
		//deflection
		Eigen::Vector3d normal = (m_A * (intersection - m_pos)).normalized();
		particle->setVelocity(particle->getVelocity() - 2 * normal*normal.dot(particle->getVelocity()));
	}
	else {
		//transmission
		if (contains(oldPos)) {
			//we are leaving this compartment and going to parent
			particle->changeCompartment(m_parentComp, timeStep);
		}
		else {
			//we are entering this compartment
			particle->changeCompartment(std::make_shared<Ellipsoid>(*this), timeStep);
		}
	}
	return false;
}




