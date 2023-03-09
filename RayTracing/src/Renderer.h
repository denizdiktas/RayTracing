#pragma once

#include "Walnut/Image.h"

#include "Camera.h"
#include "Ray.h"
#include "Scene.h"

#include <memory>
#include <glm/glm.hpp>

#include <atomic>
#include "TileBeam.h"


// PREPROCESSOR SWITCHES
#define MT
#define THREAD_LOCAL_RANDOM
#define USE_CACHED_RANDOM_NORMALS

// IMPORTANT: you have to select one of the task granularity levels below
//#define MT_TASK_GRANULARITY_PIXEL
//#define MT_TASK_GRANULARITY_ROW
#define MT_TASK_GRANULARITY_COL
//#define MT_TASK_GRANULARITY_TILE

// the following can be used only for tile-based rendering
#define USE_TILE_BEAM_INTERSECTION_TEST


class Renderer
{
public:
	struct Settings
	{
		bool Accumulate = true;
	};
public:
	Renderer();

	void OnResize(uint32_t width, uint32_t height);
	void Render(const Scene& scene, const Camera& camera);
	void CalcImageData(int x, int y);
	void UpdateImageData(int x, int y, const glm::vec4& color); // use this to set the color directly in the beam-intersection test!

	std::shared_ptr<Walnut::Image> GetFinalImage() const { return m_FinalImage; }

	void ResetFrameIndex() { m_FrameIndex = 1; }
	Settings& GetSettings() { return m_Settings; }
private:
	struct HitPayload
	{
		float HitDistance;
		glm::vec3 WorldPosition;
		glm::vec3 WorldNormal;

		int ObjectIndex;
	};

	glm::vec4 PerPixel(uint32_t x, uint32_t y); // RayGen

	HitPayload TraceRay(const Ray& ray);
	HitPayload ClosestHit(const Ray& ray, float hitDistance, int objectIndex);
	HitPayload Miss(const Ray& ray);
private:
	std::shared_ptr<Walnut::Image> m_FinalImage;
	Settings m_Settings;

	std::vector<uint32_t> m_ImageHorizontalIter, m_ImageVerticalIter;

	const Scene* m_ActiveScene = nullptr;
	const Camera* m_ActiveCamera = nullptr;

	uint32_t* m_ImageData = nullptr;
	glm::vec4* m_AccumulationData = nullptr;

	uint32_t m_FrameIndex = 1;

protected:
	const int				m_NumRandomNormals = 1024 * 1024;
	std::vector<glm::vec3>	m_RandomNormals;

	std::atomic<int>	m_GlobalThreadCount = 0; // keeps track of the total number of threads in the thread-pool
	std::vector<float>	m_TotalFrameTimePerThread; // EACH ENTRY keeps track of the TOTAL FRAME TÝME FOR EACH THREAD

	const int m_TileSizeX = 8;
	const int m_TileSizeY = 8;
	int m_NumTilesX, m_NumTilesY;
	std::vector<int>		m_TileIterX, m_TileIterY;
	std::vector<TileBeam>	m_TileBeams;
	bool m_UpdateTileBeams = true;
};