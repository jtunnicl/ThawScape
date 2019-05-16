#pragma once

#include <vector>
#include <random>
#include <numeric>
#include <algorithm>
#include "Array2D.hpp"

#define NR_END 1
#define FREE_ARG char*

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
//#define M 7
#define NSTACK 100000

#define sqrt2 1.414213562373f
#define oneoversqrt2 0.707106781186f
#define degrad 0.01745329251994330   // Convert degrees to radians; e.g. 180 * degrad = 3.14159..
#define PI 3.14159265358979
#define HALFPI = PI/2
#define fillincrement 0.01f


class time_fcn {
public:
	int year;
	int day;           // Everything works on 365 Julian Day system for now
	int hour;          // 24 hr clock
	int minute;
	int LocalTime;     // Decimal Hour
	int UT;

	int end_year;      // Model end year
};

class solar_geom {
public:
	float lattitude;
	float longitude;
	float stdmed;          // LSTM = (UTC - 7H * 15 deg)
	float declination;     // Declination of sun from equatorial plane
	float altitude;        // Sun altitude in the sky
	float azimuth;         // Compass angle of sun
	float incidence;       // Angle of sun's incidence
	float SHA;             // Solar Hour Angle is 0° at solar noon.  Since the Earth rotates 15° per hour,
		            // each hour away from solar noon corresponds to an angular motion of the sun in the sky of 15°.
		            // In the morning the hour angle is negative, in the afternoon the hour angle is positive.
};

class StreamPower
{
public:

	int lattice_size_x, lattice_size_y, duration, printinterval, printstep;
	float U, K, D, melt, timestep, ann_timestep, deltax, deltax2, thresh, thresh_diag, thresholdarea;
	float init_exposure_age, init_sed_track, init_veg;

	// new vars
	float xllcorner, yllcorner, nodata;
	std::vector<int> iup, idown, jup, jdown, topovecind, sed_vecind;
	std::vector<float> ax, ay, bx, by, cx, cy, ux, uy, rx, ry, topovec, sed_vec;
	std::vector<std::vector<float>> topo, topoold, topo2, slope, aspect, flow, flow1, flow2, flow3, 
		flow4, flow5, flow6, flow7, flow8, FA, veg, veg_old, Sed_Track, ExposureAge, ExposureAge_old;
	std::vector<std::vector<float>> solar_raster, shade_raster, I_D;
	std::vector<std::vector<std::vector<float>>> Ip_D8;     // Map of incoming solar flux, 8 directions
	Array2D<float> elevation;

	time_fcn ct;             // Current model time
	solar_geom r;

	static std::vector<float> Vector(int nl, int nh);
	static std::vector<int> IVector(int nl, int nh);
	static std::vector<std::vector<float>> Matrix(int nrl, int nrh, int ncl, int nch);			
	static std::vector<std::vector<int>> IMatrix(int nrl, int nrh, int ncl, int nch);	

	static float Ran3(std::default_random_engine& generator, std::uniform_real_distribution<float>& distribution);
	static float Gasdev(std::default_random_engine& generator, std::normal_distribution<float>& distribution);

	static void Indexx(int n, float* arr, int* indx);	// interface from old to new implementation
	static std::vector<int> Indexx(std::vector<float>& arr);	// new implementation

	static void Tridag(float a[], float b[], float c[], float r[], float u[], unsigned long n); // interface from old to new implementation
	static void Tridag(std::vector<float>& a, std::vector<float>& b, std::vector<float>& c, std::vector<float>& r, std::vector<float>& u, int n); // new implementation

	StreamPower(int nx, int ny);
	~StreamPower();

	std::vector<std::vector<float>> CreateRandomField();
	std::vector<std::vector<float>> ReadArcInfoASCIIGrid(char* fname);
	std::vector<std::vector<float>> GetTopo();


	void SetupGridNeighbors();
	void SetTopo(std::vector<std::vector<float>> t);
	void SetFA(std::vector<std::vector<float>> f);     // Set flow accumulation raster
	void Flood(); // Barnes pit filling
	void MFDFlowRoute(int i, int j); //new implementation
	void InitDiffusion();
	void HillSlopeDiffusion();
	void Avalanche(int i, int j);
	void SlopeAspect(int i, int j);
	void SunPosition();
	void SolarInflux();
	void MeltExposedIce();

	void Init(); // using new vars
	void Start();
	void PrintState(char* fname);
};

template <typename T> std::vector<T> ArrayToVector(T* a, int size)
{
	std::vector<T> v = std::vector<T>(size);
	for (int i = 0; i < size; i++)
	{
		v[i] = a[i];
	}
	return v;
}

template <typename T> std::vector<T> ArrayToVector(T* a, int size, bool fortranIndexing)
{
	std::vector<T> v;
	if (fortranIndexing)
	{
		v = std::vector<T>(size + 1);
	}
	else 
	{
		v = std::vector<T>(size);
	}

	for (int i = 0; i < size; i++)
	{
		v[i] = a[i];
	}
	return v;
}

template <typename T> void VectorToArray(std::vector<T>& v, T* a)
{
	for (int i = 0; i < v.size(); i++)
	{
		a[i] = v[i];
	}
}

// http://stackover_flow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template <typename T> std::vector<int> SortIndices(const std::vector<T>& v)
{

	// initialize original index locations
	std::vector<int> idx(v.size());
	std::iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	std::sort(idx.begin(), idx.end(), [&v](int i1, int i2) {return v[i1] < v[i2]; });

	return idx;
}

template <typename T> std::vector<int> SortFortranIndices(const std::vector<T>& v)
{

	// initialize original index locations
	std::vector<int> idx(v.size());
	std::iota(idx.begin()+1, idx.end(), 1);

	// sort indexes based on comparing values in v
	std::sort(idx.begin()+1, idx.end(), [&v](int i1, int i2) {return v[i1] < v[i2]; });

	return idx;
}
