#ifndef _STREAMPOWER_H_
#define _STREAMPOWER_H_

#include <vector>
#include <random>
#include <numeric>
#include <algorithm>
#include "global_defs.h"
#include "model_time.h"
#include "raster.h"
#include "mfd_flow_router.h"
#include "grid_neighbours.h"
#include "parameters.h"
#include "hillslope_diffusion.h"
#include "dem.h"
#include "radiation_model.h"
#include "avalanche.h"
#include "flood.h"

#define NR_END 1
#define FREE_ARG char*

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
//#define M 7
#define NSTACK 100000

#define HALFPI = PI/2
#define fillincrement 0.01


class StreamPower
{
public:

	int lattice_size_x, lattice_size_y, printstep;
	real_type deltax, deltax2;

	// new vars
    Parameters params;
	real_type xllcorner, yllcorner, nodata;
	std::vector<int> iup, idown, jup, jdown;
    DEM topo;
    Raster flow;
	Raster veg, veg_old, Sed_Track, ExposureAge, ExposureAge_old;
    MFDFlowRouter mfd_flow_router;
    GridNeighbours nebs;
    HillSlopeDiffusion hillslope_diffusion;
    RadiationModel radiation_model;
    Avalanche avalanche;
    Flood flood;

	ModelTime ct;             ///< Current model time

	static std::vector<real_type> Vector(int nl, int nh);
	static std::vector<int> IVector(int nl, int nh);
	static std::vector<std::vector<real_type>> Matrix(int nrl, int nrh, int ncl, int nch);			
	static std::vector<std::vector<int>> IMatrix(int nrl, int nrh, int ncl, int nch);	

    bool fix_random_seed;
	static real_type Ran3(std::default_random_engine& generator, std::uniform_real_distribution<real_type>& distribution);
	static real_type Gasdev(std::default_random_engine& generator, std::normal_distribution<real_type>& distribution);

    real_type channel_erosion();
    void uplift();

	StreamPower(int nx, int ny);
	~StreamPower();

	std::vector<std::vector<real_type>> CreateRandomField();


	void SetTopo();
	void SetFA();
	void InitDiffusion();

	void Init(std::string parameter_file); // using new vars
	void Start();

    std::string topo_file, fa_file, sed_file;
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

#endif
