#define _USE_MATH_DEFINES
#include <cmath>

#include <libSystem_position.h>
#include <libSystem_propagation.h>

//
// function absolute_angle
//

radians absolute_angle (const Point &positionA, const Point &positionB)
{
    // Computes the angle in radians (-PI, +PI) between two
    // positions as seen from positionA w.r.t the horizontal reference (0 radians)                  

    double xDiff = positionB.x - positionA.x;
    double yDiff = positionB.y - positionA.y;

    radians angle = atan2<radians>(yDiff, xDiff);

	return angle;
}

//
// class PositionList
//

PositionList::PositionList(const PointList &points, Site * const site, const FadeParams &fadeParams) : vector<Position*>(0)
{
    for (size_t i = 0; i < points.size(); i++)
        push_back(new Position(points[i], site, fadeParams));
}

void PositionList::clear_and_delete()
{
    while (!empty())
    {
        const Position *position = back();
        pop_back();
        delete position;
    }
}

void PositionList::reset_shadow_fading ()
{
    for (size_t i=0; i<size(); i++)
        this->at(i)->reset();
}

//
// class Position
//

Position::Position (const Point &point, const Site * const site, const FadeParams &fadeParams)
    : point          (point),
      site           (site),
      _shadowFading  (fadeParams.shadowFadingStdDev),
      _buildingFading(fadeParams.buildingFadingStdDev, fadeParams.buildingFadingMedian)
      {}
