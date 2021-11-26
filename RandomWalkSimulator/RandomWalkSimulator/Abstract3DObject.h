#ifndef ABSTRACT3DOBJECT_H
#define ABSTRACT3DOBJECT_H

#include <Eigen/Dense>

class Abstract3DObject
{
protected:
	Eigen::Vector3d m_pos;

	Abstract3DObject(double x, double y, double z);

public:
	void setPos(Eigen::Vector3d pos);

	Eigen::Vector3d getPos() { return m_pos; }
};

#endif
