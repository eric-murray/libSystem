#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include <algorithm>
#include <vector>

#include <libEDM_propagation.h>

#include <libSystem_interference.h>
#include <libSystem_propagation.h>
#include <libSystem_sites.h>

using std::min;
using std::pair;
using std::random_shuffle;
using std::vector;

class Network : public SiteList<Site> {
public:
	NetworkPropagation &propagation;

    const bool randomiseActiveNodeBs;

	Network (NetworkPropagation &propagation, const bool randomiseActiveNodeBs) : propagation(propagation), randomiseActiveNodeBs(randomiseActiveNodeBs), _numNodeBs(0), _xHi(0.0), _xLo(0.0), _yHi(0.0), _yLo(0.0) {Interference::initialise(this);}

    metres xWidth() const {return xHi() - xLo();}
    metres yWidth() const {return yHi() - yLo();}

    size_t numInterferers () const {return _interferers.size();}
    size_t numNodeBs      () const {return _numNodeBs;}
    size_t numUEs         () const;

    void add_interferer (Interferer *interferer) {_interferers.push_back(interferer);}
    void add_nodeBs     (const metres x, const metres y, const NodeBParams &nodeBParams, const radians startAzimuth = radians());
    void add_nodeBs     (Site *site, const NodeBParams &nodeBParams, const radians startAzimuth = radians());

    double compute_interference (const Station &station);
    void   create_interferers   (const size_t numInterferers, const dBm &interfererPower, const metres &height, const FadeParams &fadeParams, const bool tessellate = false, const bool generateNumInterferersPerCell = false);
    void   create_interferer    (const Position &position, const dBm &interfererPower, const metres &height, const bool tessellate);
    metres closest_interferer   (const Station &station) const;

    void randomise_active_nodeBs (const NodeB &connectedNodeB);

    // virtual functions
    virtual metres xHi() const {return _xHi;}
    virtual metres xLo() const {return _xLo;}
    virtual metres yHi() const {return _yHi;}
    virtual metres yLo() const {return _yLo;}

    virtual size_t clusterSize() const {return 1;}

    // pure virtual functions
    virtual double   deploymentArea           ()                             const = 0;
    virtual Position generate_random_position (const FadeParams &fadeParams)       = 0;

    bool check_powers       () const;
    bool connected_calls_ok () const;
    void print_sites        (string  &filename) const {ofstream file(filename.c_str()); print_sites(file); file.close();}
    void print_sites        (ostream &file = cout, const bool printTessellations = false) const;
    void print_parameters   (ostream &file = cout) const;
	void print_ues          (ostream &file = cout, const bool banner = true, const size_t leadingSpaces = 0) const;
    void reset              ();

protected:
    metres _xHi, _xLo, _yHi, _yLo;

    // pure virtual functions
    virtual PointList tessellate_points (const metres &x, const metres &y) const = 0;

private:
    ulong          _numNodeBs;
    InterfererList _interferers;
};

class CircularNetwork : public Network {
public:
    const metres radius;

    double deploymentArea() const {return M_PI * radius * radius;}

    CircularNetwork(NetworkPropagation &propagation, const metres radius, const NodeBParams &nodeBParams, const bool randomiseActiveNodeBs) : Network(propagation, randomiseActiveNodeBs), radius(radius) {add_nodeBs(metres(0.0), metres(0.0), nodeBParams);}

    Position  generate_random_position (const FadeParams &fadeParams);
    PointList tessellate_points        (const metres &x, const metres &y) const;
};

class RectangularNetwork : public Network {
public:
    metres xHi() const {return max(_xHi, _widthXCell * _numXCells);}
    metres xLo() const {return min(_xLo, metres(0.0));}
    metres yHi() const {return max(_yHi, _widthYCell * _numYCells);}
    metres yLo() const {return min(_yLo, metres(0.0));}

    double deploymentArea() const {return _numXCells * _numYCells * _widthXCell * _widthYCell;}

	RectangularNetwork(NetworkPropagation &propagation, const size_t numXCells, const size_t numYCells, const metres widthXCell, const metres widthYCell, const bool randomiseActiveNodeBs);

    void      create_interferers       (const size_t numInterferers, const dBm interfererPower, const metres height, const FadeParams &fadeParams, const bool tessellate = false);
    Position  generate_random_position (const FadeParams &fadeParams);
    void      print_parameters         (ostream &file = cout) const;
    PointList tessellate_points        (const metres &x, const metres &y) const;

private:
    typedef pair<size_t,size_t> cellid;

    const size_t _numXCells,  _numYCells;
    const metres _widthXCell, _widthYCell;

    vector<cellid> cellList, randomisedCellList;

    void reset_randomisedCellList();
};

class RegularNetwork : virtual public Network {
public:
	const metres  cellRadius;
    const size_t  numTessellations;
    const radians tessellationAngle;
    const double  tessellationDistance;
    const metres  siteToSiteDistance;

    double deploymentArea() const {return _deploymentArea;}

    RegularNetwork (NetworkPropagation &propagation, const metres siteToSiteDistance, const double cellRadius, const size_t clusterSize, const size_t numTiers, const size_t numTessellations = 1, const radians tessellationAngle = radians(), const double tessellationDistance = 0.0, const bool randomiseActiveNodeBs = false)
        : Network(propagation, randomiseActiveNodeBs), siteToSiteDistance(siteToSiteDistance), cellRadius(cellRadius), _clusterSize(clusterSize), _numTiers(numTiers), numTessellations(numTessellations), tessellationAngle(tessellationAngle), tessellationDistance(tessellationDistance), _deploymentArea(0.0) {}

    Position generate_random_position (const FadeParams &fadeParams);

    // virtual functions
    virtual size_t clusterSize() const {return _clusterSize;}

    // pure virtual functions
    virtual radians wrappedAngle    (const size_t positionId) const = 0;
    virtual double  wrappedDistance (const size_t positionId) const = 0;

protected:
    const size_t _clusterSize;
    const size_t _numTiers;

    double _deploymentArea;

    PointList tessellate_points (const metres &x, const metres &y) const;

    Site *add_site (const metres &x, const metres &y, const size_t clusterId, const FadeParams &fadeParams, const double centreOffset);

	// pure virtual functions
	virtual void  initialise (const NodeBParams &nodeBParams, const double centreOffset, const radians &startAzimuth)                                                                                                  = 0;
	virtual Site *new_site   (const metres &x, const metres &y, const size_t clusterId, const FadeParams &fadeParams, const double centreOffset) const = 0;

private:
    Point wrappedPoint (const metres x, const metres y, const size_t positionId) const {return Point(x + wrappedDistance(positionId)*cos(wrappedAngle(positionId)), y + wrappedDistance(positionId)*sin(wrappedAngle(positionId)));}
};


class RegularHexagonalNetwork : public RegularNetwork {
public:
    RegularHexagonalNetwork (NetworkPropagation &propagation, const metres siteToSiteDistance, const double cellRadius, const size_t clusterSize, const size_t numTiers, const bool randomiseActiveNodeBs = false);

protected:
	void initialise (const NodeBParams &nodeBParams, const double centreOffset, const radians &startAzimuth);

private:
    size_t  numClusters     ()                        const {return (_numTiers+1)*_numTiers*3 + 1;}
    radians wrappedAngle    (const size_t positionId) const {return tessellationAngle + PI / 3.0 * positionId;}
    double  wrappedDistance (const size_t positionId) const {return tessellationDistance;}

    void add_clusters (const metres &x, const metres &y, const NodeBParams &nodeBParams, const radians &startAzimuth, const double centreOffset);

    radians compute_tessellation_angle    (const size_t numTiers, const size_t clusterSize)                                  const;
	double  compute_tessellation_distance (const size_t numTiers, const size_t clusterSize, const double siteToSiteDistance) const;
};

class HexagonalNetwork : public RegularHexagonalNetwork {
public:
    HexagonalNetwork (NetworkPropagation &propagation, const metres siteToSiteDistance, const size_t clusterSize, const size_t numTiers, const NodeBParams &nodeBParams, const bool randomiseActiveNodeBs, const double centreOffset = 0.0, const radians startAzimuth = PI / 2.0);

private:
	HexagonalSite *new_site (const metres &x, const metres &y, const size_t clusterId, const FadeParams &fadeParams, const double centreOffset) const {return new HexagonalSite(*this, tessellate_points(x, y), clusterId, fadeParams, cellRadius, centreOffset);}
};

class CloverLeafNetwork : public RegularHexagonalNetwork {
public:
    CloverLeafNetwork (NetworkPropagation &propagation, const metres siteToSiteDistance, const size_t clusterSize, const size_t numTiers, const NodeBParams &nodeBParams, const bool randomiseActiveNodeBs, const double centreOffset = 0.0, const radians startAzimuth = PI / 2.0);

private:
	CloverLeafSite *new_site (const metres &x, const metres &y, const size_t clusterId, const FadeParams &fadeParams, const double centreOffset) const {return new CloverLeafSite(*this, tessellate_points(x, y), clusterId, fadeParams, cellRadius, centreOffset);}
};

class SquareNetwork : public RegularNetwork {
public:
    SquareNetwork (NetworkPropagation &propagation, const metres siteToSiteDistance, const size_t clusterSize, const size_t numTiers, const NodeBParams &nodeBParams, const bool randomiseActiveNodeBs, const double centreOffset = 0.0, const radians startAzimuth = radians());

private:
    size_t  numClusters     ()                        const {return (2*_numTiers + 1) * (2*_numTiers + 1);}
    radians wrappedAngle    (const size_t positionId) const {return PI * 0.25 * positionId;}
    double  wrappedDistance (const size_t positionId) const {return tessellationDistance * (positionId % 2 ? 1.0 : M_SQRT1_2);}

	void        initialise   (const NodeBParams &nodeBParams, const double centreOffset, const radians &startAzimuth);
    void        add_clusters (const metres &x, const metres &y, const NodeBParams &nodeBParams, const radians &startAzimuth);
	SquareSite *new_site     (const metres &x, const metres &y, const size_t clusterId, const FadeParams &fadeParams, const double centreOffset) const {return new SquareSite(*this, tessellate_points(x, y), clusterId, fadeParams, cellRadius, centreOffset);}
};