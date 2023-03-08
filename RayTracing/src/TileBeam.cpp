
#include "TileBeam.h"

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>

#include "Walnut/Input/Input.h"

using namespace Walnut;

void TileBeam::computeFaceNormals(const glm::vec3& d0, const glm::vec3& d1, const glm::vec3& d2, const glm::vec3& d3)
{
	faceNormals[0] = glm::normalize(glm::cross(d0, d1)); // PLANE #0: LOWER FACE 
	faceNormals[1] = glm::normalize(glm::cross(d1, d2)); // PLANE #1: RIGHT FACE 
	faceNormals[2] = glm::normalize(glm::cross(d2, d3)); // PLANE #2: UPPER FACE 
	faceNormals[3] = glm::normalize(glm::cross(d3, d0)); // PLANE #3: LEFT FACE
}

bool TileBeam::intersects(const Sphere& s, const glm::vec3& eye)
{
	const auto diff = s.Position - eye;

	// special case for the eye position and the mid-plane
	if (glm::dot(diff, -midRayDir) > 0 && glm::dot(diff, diff) > s.Radius * s.Radius)
	{
		return false;
	}

	for (const auto& n : faceNormals)
	{
		if (glm::dot(diff, n) > s.Radius)
		{
			return false;
		}
	}

	return true;
}