#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include <vector>

#include <libEDM_types.h>
#include <libEDM_library.h>
#include <libEDM_propagation.h>
#include <libEDM_random.h>

using std::vector;

// forward declarations
class FadeParams;
class Site;

//
// class Point
//
class Point {
public:
    const metres x, y;

    Point (const metres x = metres(0.0), const metres y = metres(0.0)) : x(x), y(y) {}
    Point (const Point &point                                        ) : x(point.x), y(point.y) {}

    bool  operator== (const Point &point) const {return ((point.x == x) && (point.y == y));}
    Point operator+  (const Point &point) const {return Point(x + point.x, y + point.y);}
    Point operator-  (const Point &point) const {return Point(x - point.x, y - point.y);}
    Point operator=  (const Point &point) const {return Point(point.x, point.y);}
};

//
// class PointList
//

typedef vector<Point> PointList;

//
// class Position
//

class Position {
public:
    const Point        point;
    const Site * const site;

private:
	ShadowFading _shadowFading;
    ShadowFading _buildingFading;

public:
    Position (const Point &point, const Site * const site, const FadeParams &fadeParams);

    double shadowFading  () const {return _shadowFading;}
    double buildingFading() const {return _buildingFading;}

    dB shadowFadingStdDev  () const {return _shadowFading.standardDeviation;}
    dB buildingFadingMedian() const {return _buildingFading.median;}
    dB buildingFadingStdDev() const {return _buildingFading.standardDeviation;}

    bool operator== (const Position &position) const {return (this->point.x == point.x) && (this->point.y == point.y);}

    void reset() {_shadowFading.reset(); _buildingFading.reset();}
};

//
// class PositionList
//

class PositionList : public vector<Position*> {
public:
    PositionList()                   : vector<Position*>(0) {}
    PositionList(Position &position) : vector<Position*>(0) {this->push_back(&position);}

    PositionList(const PointList &points, Site * const site, const FadeParams &fadeParams);

    const Position &position (const size_t id) const {return *at(id);}

    void clear_and_delete();
    void reset_shadow_fading();
};

//
// global functions
//

radians absolute_angle (const Point &positionA, const Point &positionB);