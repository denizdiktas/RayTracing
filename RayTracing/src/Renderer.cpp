#include "Renderer.h"

#include "Walnut/Random.h"
#include "Walnut/Timer.h"
#include "ThreadLocalRandom.h"

#include <atomic>
#include <execution>


// PREPROCESSOR SWITCHES
#define MT
#define THREAD_LOCAL_RANDOM
#define USE_CACHED_RANDOM_NORMALS

// IMPORTANT: you have to select one of the task granularity levels below
//#define MT_TASK_GRANULARITY_PIXEL
//#define MT_TASK_GRANULARITY_ROW
//#define MT_TASK_GRANULARITY_COL
#define MT_TASK_GRANULARITY_TILE
static const int tileSizeX = 8;
static const int tileSizeY = 8;

// the following can be used only for tile-based rendering
#define USE_TILE_BEAM_INTERSECTION_TEST


namespace Utils {

	static uint32_t ConvertToRGBA(const glm::vec4& color)
	{
		uint8_t r = (uint8_t)(color.r * 255.0f);
		uint8_t g = (uint8_t)(color.g * 255.0f);
		uint8_t b = (uint8_t)(color.b * 255.0f);
		uint8_t a = (uint8_t)(color.a * 255.0f);

		uint32_t result = (a << 24) | (b << 16) | (g << 8) | r;
		return result;
	}

}


namespace {
	const int numRandomNormals = 1024 * 1024;
	thread_local int currentRandomNormal = 0;
	std::vector<glm::vec3> randomNormals;
	
	std::atomic<int> globalThreadCount = 0; // keeps track of the total number of threads in the thread-pool
	std::vector<float> totalFrameTimePerThread; // EACH ENTRY keeps track of the TOTAL FRAME TÝME FOR EACH THREAD

	std::vector<int> tileIterX, tileIterY;
}

Renderer::Renderer()
{
#ifdef USE_CACHED_RANDOM_NORMALS
	randomNormals.reserve(numRandomNormals);
	for (int i = 0; i < numRandomNormals; i++)
		randomNormals.push_back(Walnut::Random::Vec3(-.5, .5));
	//randomNormals.push_back(glm::normalize(Walnut::Random::Vec3(-.5, .5)));
#endif
}

void Renderer::OnResize(uint32_t width, uint32_t height)
{
	if (m_FinalImage)
	{
		// No resize necessary
		if (m_FinalImage->GetWidth() == width && m_FinalImage->GetHeight() == height)
			return;

		m_FinalImage->Resize(width, height);
	}
	else
	{
		m_FinalImage = std::make_shared<Walnut::Image>(width, height, Walnut::ImageFormat::RGBA);
	}

	delete[] m_ImageData;
	m_ImageData = new uint32_t[width * height];

	delete[] m_AccumulationData;
	m_AccumulationData = new glm::vec4[width * height];

	m_ImageHorizontalIter.resize(width);
	m_ImageVerticalIter.resize(height);
	for (uint32_t i = 0; i < width; i++)
		m_ImageHorizontalIter[i] = i;
	for (uint32_t i = 0; i < height; i++)
		m_ImageVerticalIter[i] = i;


	// for tiled rendering
	const auto numTilesX = (width + tileSizeX - 1) / tileSizeX;
	const auto numTilesY = (height + tileSizeY - 1) / tileSizeY;
	tileIterX.resize(numTilesX);
	tileIterY.resize(numTilesY);
	for (int i = 0; i < numTilesX; i++)
		tileIterX[i] = i;
	for (int i = 0; i < numTilesY; i++)
		tileIterY[i] = i;

	// be conservative for the number of threads (max of columns and rows)
	const auto maxThreads = std::max(width, height);

	// NOTE: for MT_TASK_GRANULARITY_PIXEL, the total number of threads might seem to be much more, but it turns out not to be the case
	//       but if you suspect that it might be the case, use the following line:
	// const auto maxThreads = width * height;

	totalFrameTimePerThread.resize(maxThreads, 0);
}


#ifdef USE_TILE_BEAM_INTERSECTION_TEST
struct TileBeam
{
	glm::vec3 midRayDir;
	glm::vec3 faceNormals[4];

	void computeFaceNormals(const glm::vec3& d0, const glm::vec3& d1, const glm::vec3& d2, const glm::vec3& d3)
	{
		faceNormals[0] = glm::normalize(glm::cross(d0, d1)); // PLANE #0: LOWER FACE 
		faceNormals[1] = glm::normalize(glm::cross(d1, d2)); // PLANE #1: RIGHT FACE 
		faceNormals[2] = glm::normalize(glm::cross(d2, d3)); // PLANE #2: UPPER FACE 
		faceNormals[3] = glm::normalize(glm::cross(d3, d0)); // PLANE #3: LEFT FACE
	}

	// all beams originate at EYE, they share this point as their COMMON ORIGIN
	bool intersects(const Sphere& s, const glm::vec3& eye)
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
};
#endif

void Renderer::Render(const Scene& scene, const Camera& camera)
{
	m_ActiveScene = &scene;
	m_ActiveCamera = &camera;
	
	if (m_FrameIndex == 1)
		memset(m_AccumulationData, 0, m_FinalImage->GetWidth() * m_FinalImage->GetHeight() * sizeof(glm::vec4));


	// accTime keeps track of the total time per frame including the SCHEDULING OVERHEAD
	static float  accTime = 0.f, lastAccTime = 0.f;
	Walnut::Timer  timer;

#ifdef MT
	#ifdef MT_TASK_GRANULARITY_PIXEL
	
	std::for_each(std::execution::par, m_ImageVerticalIter.begin(), m_ImageVerticalIter.end(), [this](uint32_t y)
	{
		std::for_each(std::execution::par, m_ImageHorizontalIter.begin(), m_ImageHorizontalIter.end(), [this, y](uint32_t x)
		{
			thread_local static const int tid = globalThreadCount++;
			Walnut::Timer localTimer;

			CalcImageData(x, y);

			const auto localElapsedTime = localTimer.ElapsedMillis();
			totalFrameTimePerThread[tid] += localElapsedTime;
		});
	});
	
	#elif defined(MT_TASK_GRANULARITY_ROW)

	std::for_each(std::execution::par, m_ImageVerticalIter.begin(), m_ImageVerticalIter.end(), [this](uint32_t y)
	{
		thread_local static const int tid = globalThreadCount++;
		Walnut::Timer localTimer;
			
		for (uint32_t x = 0; x < m_FinalImage->GetWidth(); x++)
		{
			CalcImageData(x, y);
		}

		const auto localElapsedTime = localTimer.ElapsedMillis();
		totalFrameTimePerThread[tid] += localElapsedTime;
	});

	#elif defined(MT_TASK_GRANULARITY_COL)

	std::for_each(std::execution::par, m_ImageHorizontalIter.begin(), m_ImageHorizontalIter.end(), [this](uint32_t x)
	{
		thread_local static const int tid = globalThreadCount++;
		Walnut::Timer localTimer;

		for (uint32_t y = 0; y < m_FinalImage->GetHeight(); y++)
		{
			CalcImageData(x, y);
		}

		const auto localElapsedTime = localTimer.ElapsedMillis();
		totalFrameTimePerThread[tid] += localElapsedTime;
	});

	#elif defined(MT_TASK_GRANULARITY_TILE)
	
	const auto width = m_FinalImage->GetWidth();
	const auto height = m_FinalImage->GetHeight();
	//for (int ty = 0; ty < tileIterY.size(); ty++)
	std::for_each(std::execution::par, tileIterY.begin(), tileIterY.end(), [&](uint32_t ty)
	{
			std::for_each(std::execution::par, tileIterX.begin(), tileIterX.end(), [&, ty](uint32_t tx)
				{
					thread_local static const int tid = globalThreadCount++;
					Walnut::Timer localTimer;

					const auto xmin = tx * tileSizeX;
					const auto ymin = ty * tileSizeX;
					const auto xmax = std::min<uint32_t>(xmin + tileSizeX, width) - 1;
					const auto ymax = std::min<uint32_t>(ymin + tileSizeY, height) - 1;

#ifdef USE_TILE_BEAM_INTERSECTION_TEST
					// BEAM INTERSECTION TEST: check if any sphere intersects the beam, if yes then proceed with the intersection test
					//                         this could be called "TILE-BASED VF-CULLING"
					// NOTE: the beam is defined as the pyramid starting at the EYE and SPANNED by the 4 CORNER-DIRECTIONS

					// SIMPLE COMPUTATIONS FOR NOW:
					// NOTE: there is a trade-off between the check for the computation of the intersection test
					// CAUTION: all plane-normals POINT OUTWARD !!!

					// get the corner directions:
					const auto width = m_FinalImage->GetWidth();
					auto d0 = m_ActiveCamera->GetRayDirections(xmin, ymin);
					auto d1 = m_ActiveCamera->GetRayDirections(xmax, ymin);
					auto d2 = m_ActiveCamera->GetRayDirections(xmax, ymax);
					auto d3 = m_ActiveCamera->GetRayDirections(xmin, ymax);
					const auto xmid = (xmin + xmax) / 2;
					const auto ymid = (ymin + ymax) / 2;


					TileBeam tb;
					tb.computeFaceNormals(d0, d1, d2, d3);
					tb.midRayDir = m_ActiveCamera->GetRayDirections(xmid, ymid);



					// check until a sphere is intersected
					const auto& eye = m_ActiveCamera->GetPosition();
					bool intersectsAnySphere = false;
					for (auto& currentSphere : m_ActiveScene->Spheres)
					{
						bool intersectsCurrentSphere = true;
						if (tb.intersects(currentSphere, eye))
						{
							intersectsAnySphere = true;
							break;
						}
					}

					// if no objects are hit, just set all of the pixels inside the current tile to the sky-color!
					if (!intersectsAnySphere)
					{
						static const glm::vec4 red(1, 0, 0, 1);
						static const glm::vec4 skyColor = glm::vec4(0.6f, 0.7f, 0.9f, 1.0f);
						for (int y = ymin; y <= ymax; y++)
						{
							for (int x = xmin; x <= xmax; x++)
								UpdateImageData(x, y, red); // skyColor);
				}
				return;
			}
#endif

			for (int y = ymin; y <= ymax; y++)
			{
				for (int x = xmin; x <= xmax; x++)
					CalcImageData(x, y);
			}

			const auto localElapsedTime = localTimer.ElapsedMillis();
			totalFrameTimePerThread[tid] += localElapsedTime;
		});
	});
	//}

	#endif

#else

	for (uint32_t y = 0; y < m_FinalImage->GetHeight(); y++)
	{
		for (uint32_t x = 0; x < m_FinalImage->GetWidth(); x++)
		{
			CalcImageData(x, y);
		}
	}
#endif

	auto elapsed = timer.ElapsedMillis();


	m_FinalImage->SetData(m_ImageData);

	if (m_Settings.Accumulate)
		m_FrameIndex++;
	else
		m_FrameIndex = 1;


	accTime += elapsed;
	if ((accTime - lastAccTime) > 1000)
	{
		std::cout << "---------------------------------\n";
		// show all average time in each thread:
		float minAvgTime = std::numeric_limits<float>::max();
		float maxAvgTime = std::numeric_limits<float>::min();
		for (int i = 0; i < globalThreadCount; i++)
		{
			const auto currentAvgTime = totalFrameTimePerThread[i] / m_FrameIndex;
			minAvgTime = std::min(minAvgTime, currentAvgTime);
			maxAvgTime = std::max(maxAvgTime, currentAvgTime);
			std::cout << i << ": " << currentAvgTime << std::endl;
		}
		std::cout << "*** min thread avg-time = " << minAvgTime << std::endl;
		std::cout << "*** max thread avg-time = " << maxAvgTime << std::endl;
		std::cout << "*** measured avg time including overhead = " << accTime / m_FrameIndex << std::endl;
		lastAccTime = accTime;
	}
}


void Renderer::CalcImageData(int x, int y)
{
	glm::vec4 color = PerPixel(x, y);
	UpdateImageData(x,y, color);
}

void Renderer::UpdateImageData(int x, int y, const glm::vec4& color)
{
	const auto index = x + y * m_FinalImage->GetWidth();
	m_AccumulationData[index] += color;

	glm::vec4 accumulatedColor = m_AccumulationData[index];
	accumulatedColor /= (float)m_FrameIndex;

	accumulatedColor = glm::clamp(accumulatedColor, glm::vec4(0.0f), glm::vec4(1.0f));
	m_ImageData[index] = Utils::ConvertToRGBA(accumulatedColor);
}


glm::vec4 Renderer::PerPixel(uint32_t x, uint32_t y)
{
	Ray ray;
	ray.Origin = m_ActiveCamera->GetPosition();
	ray.Direction = m_ActiveCamera->GetRayDirections(x,y);
	
	glm::vec3 color(0.0f);
	float multiplier = 1.0f;

	int bounces = 5;
	for (int i = 0; i < bounces; i++)
	{
		Renderer::HitPayload payload = TraceRay(ray);
		if (payload.HitDistance < 0.0f)
		{
			glm::vec3 skyColor = glm::vec3(0.6f, 0.7f, 0.9f);
			color += skyColor * multiplier;
			break;
		}

		glm::vec3 lightDir = glm::normalize(glm::vec3(-1, -1, -1));
		float lightIntensity = glm::max(glm::dot(payload.WorldNormal, -lightDir), 0.0f); // == cos(angle)

		const Sphere& sphere = m_ActiveScene->Spheres[payload.ObjectIndex];
		const Material& material = m_ActiveScene->Materials[sphere.MaterialIndex];

		glm::vec3 sphereColor = material.Albedo;
		sphereColor *= lightIntensity;
		color += sphereColor * multiplier;

		multiplier *= 0.5f;

		ray.Origin = payload.WorldPosition + payload.WorldNormal * 0.0001f;
		ray.Direction = glm::reflect(ray.Direction,
#ifdef USE_CACHED_RANDOM_NORMALS
			payload.WorldNormal + material.Roughness * randomNormals[currentRandomNormal++ % numRandomNormals]);
#else
	#ifdef THREAD_LOCAL_RANDOM
			payload.WorldNormal + material.Roughness * Walnut::ThreadLocal::Random::Vec3(-0.5f, 0.5f));
	#else
			payload.WorldNormal + material.Roughness * Walnut::Random::Vec3(-0.5f, 0.5f));
	#endif
#endif

	}

	return glm::vec4(color, 1.0f);
}

Renderer::HitPayload Renderer::TraceRay(const Ray& ray)
{
	// (bx^2 + by^2)t^2 + (2(axbx + ayby))t + (ax^2 + ay^2 - r^2) = 0
	// where
	// a = ray origin
	// b = ray direction
	// r = radius
	// t = hit distance

	int closestSphere = -1;
	float hitDistance = std::numeric_limits<float>::max();
	for (size_t i = 0; i < m_ActiveScene->Spheres.size(); i++)
	{
		const Sphere& sphere = m_ActiveScene->Spheres[i];
		glm::vec3 origin = ray.Origin - sphere.Position;

		float a = glm::dot(ray.Direction, ray.Direction);
		float b = 2.0f * glm::dot(origin, ray.Direction);
		float c = glm::dot(origin, origin) - sphere.Radius * sphere.Radius;

		// Quadratic forumula discriminant:
		// b^2 - 4ac

		float discriminant = b * b - 4.0f * a * c;
		if (discriminant < 0.0f)
			continue;

		// Quadratic formula:
		// (-b +- sqrt(discriminant)) / 2a

		// float t0 = (-b + glm::sqrt(discriminant)) / (2.0f * a); // Second hit distance (currently unused)
		float closestT = (-b - glm::sqrt(discriminant)) / (2.0f * a);
		if (closestT > 0.0f && closestT < hitDistance)
		{
			hitDistance = closestT;
			closestSphere = (int)i;
		}
	}

	if (closestSphere < 0)
		return Miss(ray);

	return ClosestHit(ray, hitDistance, closestSphere);
}

Renderer::HitPayload Renderer::ClosestHit(const Ray& ray, float hitDistance, int objectIndex)
{
	Renderer::HitPayload payload;
	payload.HitDistance = hitDistance;
	payload.ObjectIndex = objectIndex;

	const Sphere& closestSphere = m_ActiveScene->Spheres[objectIndex];

	glm::vec3 origin = ray.Origin - closestSphere.Position;
	payload.WorldPosition = origin + ray.Direction * hitDistance;
	payload.WorldNormal = glm::normalize(payload.WorldPosition);

	payload.WorldPosition += closestSphere.Position;

	return payload;
}

Renderer::HitPayload Renderer::Miss(const Ray& ray)
{
	Renderer::HitPayload payload;
	payload.HitDistance = -1.0f;
	return payload;
}
