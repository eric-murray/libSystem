#include <libSystem_interference.h>
#include <libSystem_networks.h>
#include <libSystem_sites.h>

//
// class RegularSiteShape
//

PointList RegularSiteShape::compute_boundary_points (const double offset)
{
    PointList boundaryPoints;;

    for (size_t i=0; i<_numCorners; i++)
    {
        const metres boundaryX = centre.x + (radius * cos(_angleOfFirstCorner + i * _sectorAngle));
        const metres boundaryY = centre.y + (radius * sin(_angleOfFirstCorner + i * _sectorAngle));
        Point boundaryPoint = Point(boundaryX, boundaryY);

        boundaryPoints.push_back(boundaryPoint);
    }

    // copy first point to end, so that boundary closes
    if ( !boundaryPoints.empty() )
        boundaryPoints.push_back(boundaryPoints.front());

    return boundaryPoints;
}

metres RegularSiteShape::distanceToBoundary (const radians &angle) const
{
    radians temp = angle - _angleOfFirstCorner;

    while (temp > _sectorAngle)
	    temp -= _sectorAngle;

    while (temp < PI * 0.0)
	    temp += _sectorAngle;

    temp = (_sectorAngle/2.0 - temp);

    metres distance = radius * _distanceToBoundaryRadiusScale / cos(temp);

    return distance;
}

//
// class HexagonalSiteShape
//

Point HexagonalSiteShape::generate_random_point (const size_t sectorId) const
{
    radians randomAngle;
    metres  randomDistance = metres (0.0);
    do {
        if ( sectorId == uMAX )
            // generate position anywhere in site
            randomAngle = globalRandom.angle();
        else
            // TODO
            randomAngle = globalRandom.angle();
        randomDistance = radius * sqrt(globalRandom.sample());
    }
    while ( randomDistance > distanceToBoundary(randomAngle) );

    return Point(centre.x + randomDistance * cos(randomAngle), centre.y + randomDistance * sin(randomAngle));
}

//
// class SquareSiteShape
//

Point SquareSiteShape::generate_random_point (const size_t sectorId) const
{
    // TODO : Add case when sectorId != uMAX
    const double xOffset = (globalRandom.sample() * M_SQRT2 - M_SQRT1_2) * radius;
    const double yOffset = (globalRandom.sample() * M_SQRT2 - M_SQRT1_2) * radius;

    return Point(centre.x + xOffset, centre.y + yOffset);
}

//
// class TriHexagonalSiteShape
//

PointList CloverLeafSiteShape::compute_boundary_points(const double offset)
{
    // initialise hexagons
    _hexagons.push_back(new HexagonalSiteShape(Point(centre.x                            , centre.y + 0.50 * radius), radius / 2.0, offset, oneDegree *  90.0));
    _hexagons.push_back(new HexagonalSiteShape(Point(centre.x - 0.25 * sqrt(3.0) * radius, centre.y - 0.25 * radius), radius / 2.0, offset, oneDegree * 210.0));
    _hexagons.push_back(new HexagonalSiteShape(Point(centre.x + 0.25 * sqrt(3.0) * radius, centre.y - 0.25 * radius), radius / 2.0, offset, oneDegree * 330.0));

    PointList boundaryPoints;
    for (size_t i=0; i<_hexagons.size(); i++)
	{
        boundaryPoints.insert(boundaryPoints.end(), _hexagons[i]->boundaryPoints.begin(), _hexagons[i]->boundaryPoints.end());
		_centres.push_back(_hexagons[i]->offsetCentre);
	}

    return boundaryPoints;
}

Point CloverLeafSiteShape::generate_random_point (const size_t sectorId) const
{
    // choose a random hexagon
    size_t hexagonId;
    if ( sectorId == uMAX )
        hexagonId = globalRandom.uniform(0,2);
    else
        hexagonId = sectorId;

    return _hexagons[hexagonId]->generate_random_point();
}

//
// class Site
//

Site::Site (const Network &network, const PointList &points, const size_t clusterId, const FadeParams &fadeParams)
    : SiteShape(points[0]), id(network.numSites()), clusterId(clusterId), _network(network), positions(PositionList(points, this, fadeParams))
{}

void Site::add_sectors (ulong &numNodeBs, const NodeBParams &nodeBParams, const radians startAzimuth)
{
    radians azimuthIncrement = PI * (2.0 / nodeBParams.numSectors);
	for (size_t sectorId = 0; sectorId < nodeBParams.numSectors; sectorId++)
	{
		radians azimuth = startAzimuth + sectorId * azimuthIncrement;
        add_nodeB(numNodeBs, nodeBParams, sectorId, azimuth);
        numNodeBs++;
	}
}

void Site::add_nodeB (const ulong nodeBId, const NodeBParams &nodeBParams, const size_t sectorId, const radians azimuth)
{
    // Default case : Site uses a common carrier
	size_t carrierId = clusterId;

    // If orthogonal sectors specified, then assign different carrier to each sector, whilst maintaining cluster pattern
	if ( nodeBParams.orthogonalSectorCarriers )
		carrierId = nodeBParams.numSectors * clusterId + sectorId;

    // Exception case : If rotating sector carriers, then assume same set of carriers is used at each cluster
    // Only really makes sense for 3 sector deployments with a cluster size of 3 and orthogonal sector carriers
    if ( nodeBParams.rotateSectorCarriers && nodeBParams.orthogonalSectorCarriers && _network.clusterSize() == 3 && nodeBParams.numSectors == 3)
        carrierId = (sectorId + clusterId) % nodeBParams.numSectors;

    // create NodeB and add to list
    NodeB *nodeB = new NodeB(nodeBId, sectorId, positions, nodeBParams, carrierId, azimuth);
    push_back(nodeB);
}

NodeB *Site::nodeB (const size_t nodeBId) const
{
    return at(nodeBId);
}

size_t Site::numNodeBs() const
{
    return size();
}

void Site::randomise_active_nodeBs (const NodeB &connectedNodeB)
{
    bool connectedNodeBFound = false;

    // reset active flag for all sectors
    for (size_t i=0; i < numNodeBs(); i++)
    {
        nodeB(i)->isActive = false;
        if ( nodeB(i)->id == connectedNodeB.id )
        {
            connectedNodeBFound = true;
            nodeB(i)->isActive = true;
        }
    }

    if ( !connectedNodeBFound )
        // set random nodeB as active
        nodeB(globalRandom.uniform(0, static_cast<int>(numNodeBs())-1))->isActive = true;
}

void Site::reset ()
{
    positions.reset_shadow_fading();
    NodeBList::reset();
}