#include "Renderer.h"

#include "Walnut/Random.h"
#include "Walnut/Timer.h"
#include "ThreadLocalRandom.h"

#include <execution>



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


Renderer::Renderer()
{
#ifdef USE_CACHED_RANDOM_NORMALS
	m_RandomNormals.reserve(m_NumRandomNormals);
	for (int i = 0; i < m_NumRandomNormals; i++)
		m_RandomNormals.push_back(Walnut::Random::Vec3(-.5, .5));
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
	m_NumTilesX = (width + m_TileSizeX - 1) / m_TileSizeX;
	m_NumTilesY = (height + m_TileSizeY - 1) / m_TileSizeY;
	m_TileIterX.resize(m_NumTilesX);
	m_TileIterY.resize(m_NumTilesY);
	for (int i = 0; i < m_NumTilesX; i++)
		m_TileIterX[i] = i;
	for (int i = 0; i < m_NumTilesY; i++)
		m_TileIterY[i] = i;
	m_UpdateTileBeams = true;

	// be conservative for the number of threads (max of columns and rows)
	const auto maxThreads = std::max(width, height);

	// NOTE: for MT_TASK_GRANULARITY_PIXEL, the total number of threads might seem to be much more, but it turns out not to be the case
	//       but if you suspect that it might be the case, use the following line:
	// const auto maxThreads = width * height;

	m_TotalFrameTimePerThread.resize(maxThreads, 0);
}




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
	mtTaskGranularityPixel();	

	#elif defined(MT_TASK_GRANULARITY_ROW)
	mtTaskGranularityRow();

	#elif defined(MT_TASK_GRANULARITY_COL)
	mtTaskGranularityCol();

	#elif defined(MT_TASK_GRANULARITY_TILE)
	mtTaskGranularityTile();

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
		for (int i = 0; i < m_GlobalThreadCount; i++)
		{
			const auto currentAvgTime = m_TotalFrameTimePerThread[i] / m_FrameIndex;
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
void Renderer::mtTaskGranularityPixel()
{
	std::for_each(std::execution::par, m_ImageVerticalIter.begin(), m_ImageVerticalIter.end(), [this](uint32_t y)
		{
			std::for_each(std::execution::par, m_ImageHorizontalIter.begin(), m_ImageHorizontalIter.end(), [this, y](uint32_t x)
				{
					thread_local static const int tid = m_GlobalThreadCount++;
					Walnut::Timer localTimer;

					CalcImageData(x, y);

					const auto localElapsedTime = localTimer.ElapsedMillis();
					m_TotalFrameTimePerThread[tid] += localElapsedTime;
				}
			);
		}
	);
}
void Renderer::mtTaskGranularityRow()
{
	std::for_each(std::execution::par, m_ImageVerticalIter.begin(), m_ImageVerticalIter.end(), [this](uint32_t y)
		{
			thread_local static const int tid = m_GlobalThreadCount++;
			Walnut::Timer localTimer;

			for (uint32_t x = 0; x < m_FinalImage->GetWidth(); x++)
			{
				CalcImageData(x, y);
			}

			const auto localElapsedTime = localTimer.ElapsedMillis();
			m_TotalFrameTimePerThread[tid] += localElapsedTime;
		}
	);
}
void Renderer::mtTaskGranularityCol()
{
	std::for_each(std::execution::par, m_ImageHorizontalIter.begin(), m_ImageHorizontalIter.end(), [this](uint32_t x)
		{
			thread_local static const int tid = m_GlobalThreadCount++;
			Walnut::Timer localTimer;

			for (uint32_t y = 0; y < m_FinalImage->GetHeight(); y++)
			{
				CalcImageData(x, y);
			}

			const auto localElapsedTime = localTimer.ElapsedMillis();
			m_TotalFrameTimePerThread[tid] += localElapsedTime;
		}
	);
}
void Renderer::mtTaskGranularityTile()
{
	const auto width = m_FinalImage->GetWidth();
	const auto height = m_FinalImage->GetHeight();

#ifdef USE_TILE_BEAM_INTERSECTION_TEST
	// RECOMPUTE ALL TILE-BEAMS
	if (m_UpdateTileBeams)
	{
		m_TileBeams.resize(m_NumTilesX * m_NumTilesY);
		for (int ty = 0; ty < m_NumTilesY; ty++)
		{
			for (int tx = 0; tx < m_NumTilesX; tx++)
			{
				auto& tb = m_TileBeams[tx + ty * m_NumTilesX];

				const auto xmin = tx * m_TileSizeX;
				const auto ymin = ty * m_TileSizeY;
				const auto xmax = std::min<uint32_t>(xmin + m_TileSizeX, width) - 1;
				const auto ymax = std::min<uint32_t>(ymin + m_TileSizeY, height) - 1;

				// get the corner directions:
				const auto width = m_FinalImage->GetWidth();
				auto d0 = m_ActiveCamera->GetRayDirections(xmin, ymin);
				auto d1 = m_ActiveCamera->GetRayDirections(xmax, ymin);
				auto d2 = m_ActiveCamera->GetRayDirections(xmax, ymax);
				auto d3 = m_ActiveCamera->GetRayDirections(xmin, ymax);
				const auto xmid = (xmin + xmax) / 2;
				const auto ymid = (ymin + ymax) / 2;

				tb.computeFaceNormals(d0, d1, d2, d3);
				tb.midRayDir = m_ActiveCamera->GetRayDirections(xmid, ymid);
}
		}
		m_UpdateTileBeams = false;
	}
#endif


	//for (int ty = 0; ty < tileIterY.size(); ty++)
	std::for_each(std::execution::par, m_TileIterY.begin(), m_TileIterY.end(), [&](uint32_t ty)
		{
			std::for_each(std::execution::par, m_TileIterX.begin(), m_TileIterX.end(), [&, ty](uint32_t tx)
				{
					thread_local static const int tid = m_GlobalThreadCount++;
					Walnut::Timer localTimer;

					const auto xmin = tx * m_TileSizeX;
					const auto ymin = ty * m_TileSizeY;
					const auto xmax = std::min<uint32_t>(xmin + m_TileSizeX, width) - 1;
					const auto ymax = std::min<uint32_t>(ymin + m_TileSizeY, height) - 1;

#ifdef USE_TILE_BEAM_INTERSECTION_TEST
					// BEAM INTERSECTION TEST: check if any sphere intersects the beam, if yes then proceed with the intersection test
					//                         this could be called "TILE-BASED VF-CULLING"
					// NOTE: the beam is defined as the pyramid starting at the EYE and SPANNED by the 4 CORNER-DIRECTIONS

					// SIMPLE COMPUTATIONS FOR NOW:
					// NOTE: there is a trade-off between the check for the computation of the intersection test
					// CAUTION: all plane-normals POINT OUTWARD !!!

					auto& tb = m_TileBeams[tx + ty * m_NumTilesX];

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
					m_TotalFrameTimePerThread[tid] += localElapsedTime;
				});
		});
	//}

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

		// PERTURB the NORMAL by an amount proportonal to the material roughness
#ifdef USE_CACHED_RANDOM_NORMALS
		static thread_local int	currentRandomNormal = 0;
		auto perturbedNormal = payload.WorldNormal + material.Roughness * m_RandomNormals[currentRandomNormal++ % m_NumRandomNormals];
#else
	#ifdef THREAD_LOCAL_RANDOM
		auto perturbedNormal = payload.WorldNormal + material.Roughness * Walnut::ThreadLocal::Random::Vec3(-0.5f, 0.5f);
	#else
		auto perturbedNormal = payload.WorldNormal + material.Roughness * Walnut::Random::Vec3(-0.5f, 0.5f);
	#endif
#endif

		ray.Origin = payload.WorldPosition + payload.WorldNormal * 0.0001f;
		ray.Direction = glm::reflect(ray.Direction, perturbedNormal);
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
