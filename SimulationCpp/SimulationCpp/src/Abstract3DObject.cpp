#include "Abstract3DObject.h"

Abstract3DObject::Abstract3DObject(double x, double y, double z)
	: m_pos{ Eigen::Vector3d{x, y, z} }
{
}

void Abstract3DObject::setPos(Eigen::Vector3d pos)
{
	m_pos = pos;
}