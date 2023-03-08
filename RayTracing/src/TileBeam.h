#pragma once

#include <glm/glm.hpp>
#include "Scene.h"


struct TileBeam
{
	glm::vec3 midRayDir;
	glm::vec3 faceNormals[4];


	void computeFaceNormals(const glm::vec3& d0, const glm::vec3& d1, const glm::vec3& d2, const glm::vec3& d3);

	// all beams originate at EYE, they share this point as their COMMON ORIGIN
	bool intersects(const Sphere& s, const glm::vec3& eye);
};