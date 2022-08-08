#pragma once

#pragma warning(disable:4250)

#define _USE_MATH_DEFINES
#include <cmath>

#include <libEDM_library.h>

#include <libSystem_stations.h>

// forward declarations
class Network;
class RegularNetwork;

using std::cerr;

class SiteShape {
public:
    const Point centre;

    SiteShape (const Point &centre) : centre(centre) {}

    // virtual functions
    virtual metres xLo() const {return centre.x;}
    virtual metres xHi() const {return centre.x;}
    virtual metres yLo() const {return centre.y;}
    virtual metres yHi() const {return centre.y;}

    // pure virtual functions
    virtual double    area()                                              const {abort(); return -1.0;}
	virtual PointList centres()                                           const {abort(); return PointList(1, centre);}    
    virtual Point     generate_random_point(const size_t sectorId = uMAX) const {abort(); return centre;}
};

class RegularSiteShape : virtual public SiteShape {
public:
    const metres    radius;
    const Point     offsetCentre;
          PointList boundaryPoints;

    metres xLo() const {return centre.x + _xOffsetLo;}
    metres xHi() const {return centre.x + _xOffsetHi;}
    metres yLo() const {return centre.y + _yOffsetLo;}
    metres yHi() const {return centre.y + _yOffsetHi;}

    RegularSiteShape (const Point   &centre,
                      const metres   radius,
                      const size_t   numCorners,
                      const double   distanceToBoundaryRadiusScale,
                      const radians  angleOfFirstCorner,
                      const metres   xOffsetLo,
                      const metres   xOffsetHi,
                      const metres   yOffsetLo,
                      const metres   yOffsetHi,
					  const double   offset)
        : SiteShape                     (centre),
          radius                        (radius),
          _numCorners                   (numCorners),
          _distanceToBoundaryRadiusScale(distanceToBoundaryRadiusScale),
          _angleOfFirstCorner           (angleOfFirstCorner),
          _sectorAngle                  (PI * (2.0 / numCorners)),
          _xOffsetLo                    (xOffsetLo),
          _xOffsetHi                    (xOffsetHi),
          _yOffsetLo                    (yOffsetLo),
          _yOffsetHi                    (yOffsetHi),
		  offsetCentre                  (Point(centre.x - offset * radius * cos(angleOfFirstCorner), centre.y - offset * radius * sin(angleOfFirstCorner)))
          {}

protected:
    const size_t  _numCorners;
    const radians _sectorAngle;
    const radians _angleOfFirstCorner;
    const double  _distanceToBoundaryRadiusScale;

    metres distanceToBoundary (const radians &angle) const;

    // virtual functions
    virtual PointList compute_boundary_points(const double offset);

private:
    const metres _xOffsetLo, _xOffsetHi, _yOffsetLo, _yOffsetHi;
};

class HexagonalSiteShape : public RegularSiteShape {
public:
    double area() const {return 1.5 * sqrt(3.0) * sqr(radius);}

    HexagonalSiteShape (const Point &centre, const metres radius, const double offset, const radians angleOfFirstCorner = radians())
        : RegularSiteShape(centre, radius, 6, sqrt(3.0)/2.0, angleOfFirstCorner, -radius, radius, -radius * 0.5 * sqrt(3.0), radius * 0.5 * sqrt(3.0), offset), SiteShape(centre)
    {boundaryPoints = compute_boundary_points(offset);}

	PointList centres() const {return PointList(1, centre);}
		
	Point generate_random_point (const size_t sectorId = uMAX) const;
};

class SquareSiteShape : public RegularSiteShape {
public:
    double area() const {return 2.0 * sqr(radius);}

    SquareSiteShape (const Point &centre, const metres radius, const double offset)
        : RegularSiteShape(centre, radius, 4, M_SQRT1_2, PI / 4.0, -radius * M_SQRT1_2, radius * M_SQRT1_2, -radius * M_SQRT1_2, radius * M_SQRT1_2, 0.0), SiteShape(centre)
    {boundaryPoints = compute_boundary_points(offset);}

	PointList centres() const {return PointList(1, centre);}

protected:
	Point generate_random_point (const size_t sectorId = uMAX) const;
};

class CloverLeafSiteShape : public RegularSiteShape {
public:
    double area() const {return 1.125 * sqrt(3.0) * sqr(radius);}

    CloverLeafSiteShape (const Point &centre, const metres radius, const double offset)
        : RegularSiteShape(centre, radius, 3, 0.0, radians(), -radius * 0.5 * sqrt(3.0), radius * 0.5 * sqrt(3.0), -radius * 0.75, radius, offset), SiteShape(centre)
          {boundaryPoints = compute_boundary_points(offset);}

	PointList centres() const {return _centres;}

protected:
	Point generate_random_point (const size_t sectorId = uMAX) const;

private:
    vector<HexagonalSiteShape*> _hexagons;
	PointList                   _centres;

    PointList compute_boundary_points(const double offset);
};

class Site : public NodeBList, virtual public SiteShape {
public:
    const size_t       id;
    const size_t       clusterId;
          PositionList positions;

    Site (const Network &network, const PointList &points, const size_t clusterId, const FadeParams &fadeParams);
    ~Site() {positions.clear_and_delete();}

    NodeB  *nodeB    (const size_t nodeBId) const;
    size_t  numNodeBs()                     const;

    void add_sectors (ulong &numNodeBs, const NodeBParams &nodeBParams, const radians startAzimuth = radians());
    void add_nodeB   (const ulong nodeBId, const NodeBParams &nodeBParams, const size_t sectorId, const radians azimuth = radians());

    Position generate_random_position (const FadeParams &fadeParams)                        const {return Position(generate_random_point(),         this, fadeParams);}
    Position generate_random_position (const size_t sectorId, const FadeParams &fadeParams) const {return Position(generate_random_point(sectorId), this, fadeParams);}

    void randomise_active_nodeBs (const NodeB &connectedNodeB);

    void reset();

protected:
    const Network &_network;
};

template <class SiteShapeType>
class RegularSite : public Site, public SiteShapeType {
public:
    RegularSite (const Network &network, const PointList &points, const size_t clusterId, const FadeParams &fadeParams, const metres cellRadius, const double offset = 0.0)
        : Site(network, points, clusterId, fadeParams), SiteShapeType(centre, cellRadius, offset), SiteShape(points[0]) {}
};

typedef RegularSite<HexagonalSiteShape>  HexagonalSite;
typedef RegularSite<SquareSiteShape>     SquareSite;
typedef RegularSite<CloverLeafSiteShape> CloverLeafSite;

template <class Site_T>
class SiteList {
public:
    Site_T *site (const size_t siteId) const {return _sites.at(siteId);}

    size_t numSites() const {return _sites.size();}

    void shuffle       ()             {random_shuffle(_sites.begin(), _sites.end());}
    void push_back_site(Site_T *site) {_sites.push_back(site);}

private:
    PointerVector<Site_T> _sites;
};

typedef SiteList<HexagonalSite>  HexagonalSiteList;
typedef SiteList<SquareSite>     SquareSiteList;
typedef SiteList<CloverLeafSite> CloverLeafSiteList;