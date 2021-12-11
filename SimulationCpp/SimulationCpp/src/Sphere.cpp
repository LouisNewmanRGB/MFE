#include "Sphere.h"
#include "Simulation.h"
#include "Util.h"

#include <iostream>

Sphere::Sphere(double x, double y, double z, double T2, double diffusivity, double permeability, double radius)
	:AbstractCompartment{x, y, z, T2, diffusivity},
	m_permeability{ permeability },
	m_radius{ radius }
{
}

bool Sphere::contains(Eigen::Vector3d pos) {
	return (pos - m_pos).norm() < m_radius;
}

double Sphere::findIntersection(std::array<Eigen::Vector3d, 2> ray, double maxDistance) {
	//ray direction vector must be normalized!
	Eigen::Vector3d L = m_pos - ray[0];
	if (L.norm() - m_radius <= maxDistance) {
		double b = -2 * L.dot(ray[1]);
		double c = L.squaredNorm() - m_radius*m_radius;
		double delta = b * b - 4 * c;	
		if (delta > 0) {
			double t = (-b - sqrt(delta)) / 2;
			if (t > 0) {
				return t;
			}
			else {
				t = (-b + sqrt(delta)) / 2;
				if (t > 0) {
					return t;
				}
			}
		}
		else if (delta == 0) {
			double t = -b / 2;
			if (t > 0) {
				return t;
			}
		}
	}
	//no intersection
	return std::numeric_limits<double>::infinity();
}

bool Sphere::collide(std::shared_ptr<AbstractParticle3D> particle, Eigen::Vector3d oldPos, Eigen::Vector3d intersection, double timeStep) {
	double transmissionProba = 4 * m_permeability / particle->getSpeed();
	if (transmissionProba > 1) {
		std::cout << "Error: greater than one transmission probability!" << std::endl;
	}
	assert(transmissionProba <= 1);

	if (Util::getRandomUniform(0, 1) >= transmissionProba) {
		//deflection
		Eigen::Vector3d normal = (intersection - m_pos).normalized();
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
			particle->changeCompartment(std::make_shared<Sphere>(*this), timeStep);
		}
	}
	return false;
}
