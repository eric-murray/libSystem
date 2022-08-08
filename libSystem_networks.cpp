#define _USE_MATH_DEFINES

#include <algorithm>
#include <cmath>
#include <iostream>

#include <libSystem_interference.h>
#include <libSystem_networks.h>
#include <libSystem_propagation.h>
#include <libSystem_sites.h>

using std::atan2;
using std::cout;
using std::endl;
using std::random_shuffle;

//
// class Network
//

void Network::add_nodeBs (const metres x, const metres y, const NodeBParams &nodeBParams, const radians startAzimuth)
{
    // check if site already exists for NodeB position
    Site *site = NULL;
    for (size_t i=0; i<numSites(); i++)
        if ( x == this->site(i)->centre.x &&  y == this->site(i)->centre.y)
        {
            // site already exists
            site = this->site(i);
            break; // for
        }

    if ( site == NULL )
    {
        // iterated through entire list without finding a match, so create new site
        site = new Site(*this, tessellate_points(x, y), 0, nodeBParams.fadeParams);
        push_back_site(site);
    }

    add_nodeBs(site, nodeBParams, startAzimuth);
}

void Network::add_nodeBs (Site *site, const NodeBParams &nodeBParams, const radians startAzimuth)
{
    site->add_sectors(_numNodeBs, nodeBParams, startAzimuth);

    // update _xLo, _xHi, _yLo, _yHi
    _xLo = min(_xLo, site->xLo());
    _xHi = max(_xHi, site->xHi());
    _yLo = min(_yLo, site->yLo());
    _yHi = max(_yHi, site->yHi());
}

bool Network::check_powers() const
{
    cout << "NodeBId  TxPower  ComputedTxPower" << endl;

    for (size_t i=0; i<numSites(); i++)
        if ( !site(i)->check_powers(false) )
            return false;

    cout << endl;
    return true;
}

metres Network::closest_interferer (const Station &station) const
{
    metres closestDistance = metres::infinity();
    for (size_t i=0; i<_interferers.size(); i++)
    {
        const metres interfererDistance = station.distance(_interferers[i]->position.point);
        if (interfererDistance < closestDistance)
            closestDistance = interfererDistance;
    }

    return closestDistance;
}

double Network::compute_interference (const Station &station)
{
    double interference = 0.0;
    for (size_t i=0; i<_interferers.size(); i++)
    {
		const Interferer *interferer = _interferers[i];

        const double interfererPower = interferer->txPower();

        interference += interfererPower * propagation.pathgain(*interferer, station).pathgain;
    }
    return interference;
}

bool Network::connected_calls_ok () const
{
    for (size_t i=0; i<numSites(); i++)
        for (size_t j=0; j<site(i)->numNodeBs(); j++)
        {
            NodeB * nodeB = site(i)->nodeB(j);
            if (nodeB->numConnectedUEs() >= nodeB->maxConnectedCalls)
                return false;
        }
    return true;
}

void Network::create_interferers (const size_t numInterferers, const dBm &interfererPower, const metres &height, const FadeParams &fadeParams, const bool tessellate, const bool generateNumInterferersPerCell)
{
	for (size_t i=0; i<numSites(); i++)
    {
        Site *site = this->site(i);

		for (size_t j=0; j<site->numNodeBs(); j++)
        {
            for (size_t i=0; i < numInterferers; i++)
            {
                Position *randomPosition;
                if ( generateNumInterferersPerCell )
                {
                    const size_t sectorId = site->nodeB(j)->sectorId;
                    randomPosition = &site->generate_random_position(sectorId, fadeParams);

                    assert( j == sectorId );
                }
                else
                    randomPosition = &generate_random_position(fadeParams);

                create_interferer(*randomPosition, interfererPower, height, tessellate);
            }

            if ( !generateNumInterferersPerCell )
                return;
        }
    }
}

void Network::create_interferer (const Position &position, const dBm &interfererPower, const metres &height, const bool tessellate)
{
    PositionList positions;
    positions.push_back(new Position(position));

    if ( tessellate )
    {
        const PointList points = tessellate_points(position.point.x, position.point.y);

        for (size_t i=0; i<points.size(); i++)
            positions.push_back(new Position(points[i], position.site, FadeParams(position.shadowFadingStdDev(), position.buildingFadingMedian(), position.buildingFadingStdDev())));
    }

    Interferer *interferer = new Interferer(_interferers.size(), positions, height, interfererPower);
    _interferers.push_back(interferer);
}

void Network::print_parameters (ostream &file) const
{
    file << "Number of interferers in deployment = " << numInterferers() << endl;

    if ( numInterferers() > 0 )
    {
        const double interferencePower = _interferers.front()->txPower();
        file << "Interference power in dBm = " << watts2dBm(interferencePower) << endl;
    }
}

void Network::print_sites(ostream &file, const bool printTessellations) const
{
    file << "SiteId  " << Site::nodeBBanner.c_str() << endl;
    for (size_t i=0; i<numSites(); i++)
    {
        file << setw(6) << site(i)->id << "  ";
        site(i)->print(file, false, 8, printTessellations);
    }

    _interferers.print(file);
}

size_t Network::numUEs() const
{
    size_t numUEs = 0;

	for (size_t i=0; i<numSites(); i++)
		for (size_t j=0; j<site(i)->numNodeBs(); j++)
            numUEs += site(i)->nodeB(j)->numConnectedUEs();

    return numUEs;
}

void Network::print_ues (ostream &file, const bool banner, const size_t leadingSpaces) const
{
	bool bannerPrinted = false;
    file << fixed << setprecision(1);

	for (size_t i=0; i<numSites(); i++)
		for (size_t j=0; j<site(i)->numNodeBs(); j++)
			for (size_t k=0; k<site(i)->nodeB(j)->numConnectedUEs(); k++)
			{
			    if ( !bannerPrinted && banner )
				{
					file << "  Id       X       Y  NodeB  Pathloss  NoiseFloor  InterfererNoise  CCNoise  EffectiveNF  RcvdPwr  TxSignalPwr  RxSignalPwr    C/I  CCRcvdPwr  CCC/I  Eb/No" << endl;
					bannerPrinted = true;
				}

				const UE *ue = site(i)->nodeB(j)->ue(k);
		        ue->print(file);
			}

    if ( banner )
        file << endl;
}

void Network::randomise_active_nodeBs (const NodeB &connectedNodeB)
{
    if ( randomiseActiveNodeBs )
        for (size_t i=0; i < numSites(); i++)
            site(i)->randomise_active_nodeBs(connectedNodeB);
}

void Network::reset()
{
    // Delete all interferers
    _interferers.clearAndDelete();

    // Remove all UEs from nodeB connection lists
    for (size_t i=0; i<numSites(); i++)
        site(i)->reset();

    // reset interference calculation
    Interference::reset();
}

//
// class CircularNetwork
//

Position CircularNetwork::generate_random_position (const FadeParams &fadeParams)
{
    radians randomAngle    = globalRandom.angle();
    metres  randomDistance = radius * sqrt(globalRandom.sample());

    metres x = randomDistance * sin(randomAngle);
    metres y = randomDistance * cos(randomAngle);

    return Position(Point(x,y), NULL, fadeParams);
}

PointList CircularNetwork::tessellate_points (const metres &x, const metres &y) const
{
    // don't tessellate
    return PointList(1, Point(x, y));
}

//
// class RectangularNetwork
//

RectangularNetwork::RectangularNetwork (NetworkPropagation &propagation, const size_t numXCells, const size_t numYCells, const metres widthXCell, const metres widthYCell, const bool randomiseActiveNodeBs)
        : Network(propagation, randomiseActiveNodeBs), _numXCells(numXCells), _numYCells(numYCells), _widthXCell(widthXCell), _widthYCell(widthYCell)
{
    for (size_t x=0; x<numXCells; x++)
        for (size_t y=0; y<numYCells; y++)
            cellList.push_back(cellid(x,y));
}

void RectangularNetwork::create_interferers (const size_t numInterferers, const dBm interfererPower, const metres height, const FadeParams &fadeParams, const bool tessellate)
{
    reset_randomisedCellList();
    Network::create_interferers(numInterferers, interfererPower, height, fadeParams, tessellate);
    reset_randomisedCellList();
}

Position RectangularNetwork::generate_random_position (const FadeParams &fadeParams)
{
    if (randomisedCellList.size() == 0)
        reset_randomisedCellList();

    cellid cell = randomisedCellList.back();
    randomisedCellList.pop_back();

    // generate random position in cell for interferer
    metres x = _widthXCell * (cell.first  + globalRandom.sample());
    metres y = _widthYCell * (cell.second + globalRandom.sample());

    return Position(Point(x,y), NULL, fadeParams);
}

void RectangularNetwork::print_parameters (ostream &file) const
{
    file << "Number of cells in X dimension = "              << _numXCells                                                        << endl;
    file << "Number of cells in Y dimension = "              << _numYCells                                                        << endl;
    file << "Width of each cell in X dimension (metres) = "  << _widthXCell                                                       << endl;
    file << "Width of each cell in Y dimension (metres) = "  << _widthYCell                                                       << endl;
    file << "Fraction of cells containing an interferer = "  << static_cast<double>(numInterferers()) / (_numXCells * _numYCells) << endl;

    Network::print_parameters(file);
}

void RectangularNetwork::reset_randomisedCellList ()
{
    if ( randomisedCellList.size() != cellList.size() )
        randomisedCellList = cellList;

    random_shuffle(randomisedCellList.begin(), randomisedCellList.end());
}

PointList RectangularNetwork::tessellate_points (const metres &x, const metres &y) const
{
    // don't tessellate
    Point point = Point(x, y);
    PointList points;
    points.push_back(point);
    return points;
}

//
// class RegularNetwork
//

Site* RegularNetwork::add_site (const metres &x, const metres &y, const size_t clusterId, const FadeParams &fadeParams, const double centreOffset)
{
    Site *site = new_site(x, y, clusterId, fadeParams, centreOffset);
	push_back_site(site);

    _deploymentArea += site->area();

    return site;
}

Position RegularNetwork::generate_random_position (const FadeParams &fadeParams)
{
    // choose random site within network
    Site *site = this->site(globalRandom.uniform(0, static_cast<int>(numSites())-1));

    // get site to generate random position
    return site->generate_random_position(fadeParams);
}

PointList RegularNetwork::tessellate_points (const metres &x, const metres &y) const
{
    PointList points;
    points.push_back(Point(x, y));

    for (size_t positionId = 1; positionId <= numTessellations; positionId++)
        points.push_back(wrappedPoint(x, y, positionId));

    return points;
}

//
// class RegularHexagonalNetwork
//

RegularHexagonalNetwork::RegularHexagonalNetwork (NetworkPropagation &propagation, const metres siteToSiteDistance, const double cellRadius, const size_t clusterSize, const size_t numTiers, const bool randomiseActiveNodeBs)
    : RegularNetwork(propagation, siteToSiteDistance, cellRadius, clusterSize, numTiers, 6, compute_tessellation_angle(numTiers, clusterSize), compute_tessellation_distance(numTiers, clusterSize, siteToSiteDistance)), Network(propagation, randomiseActiveNodeBs) {}

void RegularHexagonalNetwork::add_clusters (const metres &x, const metres &y, const NodeBParams &nodeBParams, const radians &startAzimuth, const double centreOffset)
{
    Site *site = add_site(x, y, 0, nodeBParams.fadeParams, centreOffset);
    add_nodeBs(site, nodeBParams, startAzimuth);

    double  distance;
    radians angle;
    for (size_t clusterId = 1; clusterId < _clusterSize; clusterId++)
	{
        // compute distance and angle to nodeBs in cluster
	    switch (clusterId)
        {
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
            case 6:
		        distance = 1.0;
		        angle    = PI/3.0 * (clusterId-1) + PI/6.0;
                break;

            case 7:
            case 8:
            case 9:
            case 10:
            case 11:
            case 12:
		        distance = sqrt(3.0);
		        angle    = PI/3.0 * (clusterId-7);
		        if ( _clusterSize == 12 )
		            switch (clusterId)
                    {
			            case 10 :
			                distance = 2.0;
			                angle    = -PI/6.0;
                            break;
			            case 11 :
			                distance = 2.0;
			                angle    = PI/6.0;
                            break;
                        default :
			                // no adjustment required
                            break;
                    }
                break;

            case 13:
            case 14:
            case 15:
            case 16:
            case 17:
            case 18:
		        distance = 2.0;
		        angle    = PI/3.0 * (clusterId-13) + PI/6.0;
                break;
        }

        const metres clusterX = x + distance * siteToSiteDistance * cos(angle);
        const metres clusterY = y + distance * siteToSiteDistance * sin(angle);

        Site *site = add_site(clusterX, clusterY, clusterId, nodeBParams.fadeParams, centreOffset);
        add_nodeBs(site, nodeBParams, startAzimuth);
    }
}

radians RegularHexagonalNetwork::compute_tessellation_angle (const size_t numTiers, const size_t clusterSize) const
{
	switch (clusterSize)
	{
	case 1:
	case 4:
	case 9:
		return atan<radians>(1.0 / (sqrt(3.0) * (1 + 2*numTiers)));
		break;

	default:
		switch (numTiers)
		{
		case 0:
			switch (clusterSize)
			{
				case 3:
				case 12:
					return radians();       						 // 0^2 + 3*1^2 = 3 = numSites
					break;
				case 7:
					return atan<radians>(2.0 / (1.0*sqrt(3.0)));     // 2^2 + 3*1^2 = 7 = numSites
					break;
				case 13: 
					return atan<radians>(1.0 / (2.0*sqrt(3.0)));     // 1^2 + 3*2^2 = 13 = numSites
					break;
				case 19:
					return atan<radians>(4.0 / (1.0*sqrt(3.0)));     // 4^2 + 3*1^2 = 19 = numSites
					break;
			}
			break;
		case 1:
			switch (clusterSize)
			{
				case 3:
				case 12:
					return atan<radians>(3.0 / (2.0*sqrt(3.0)));     // 3^2 + 3*2^3 =  21 = numSites
					break;
				case 7:
					return PI * 0.5;                 				 // 7^2 + 3*0^2 =  49 = numSites
					break;
				case 13:
					return atan<radians>(8.0 / (3.0*sqrt(3.0)));     // 8^2 + 3*3^2 =  91 = numSites
					break;
				case 19:
					return atan<radians>(5.0 / (6.0*sqrt(3.0)));     // 5^2 + 3*6^2 = 133 = numSites
					break;
			}
			break;
		case 2:
			switch (clusterSize)
			{
				case 3:
				case 12:
					return atan<radians>( 3.0 / ( 4.0*sqrt(3.0)));   //  3^2 + 3* 4^2 =  57 = numSites
					break;
				case 7:
					return atan<radians>(13.0 / (11.0*sqrt(3.0)));   // 13^2 + 3*11^2 = 532 = numSites * 4
					break;
				case 13:
					return atan<radians>(10.0 / ( 7.0*sqrt(3.0)));   // 10^2 + 3* 7^2 = 247 = numSites
					break;
				case 19:
					return PI * 0.5;                                 // 19^2 + 3* 0^2 = 361 = numSites
					break;
			}
			break;
		case 3:
			switch (clusterSize)
			{
				case 3:
				case 12:
					return atan<radians>( 6.0 / ( 5.0*sqrt(3.0)));	 //  6^2 + 3*5^2 = 111 = numSites
					break;
				case 7:
					return atan<radians>(23.0 / (13.0*sqrt(3.0)));	 // 23^2 + 3*13^2 = 1036 = numSites * 4
					break;
				case 13:
					return atan<radians>(17.0 / ( 8.0*sqrt(3.0)));	 // 17^2 + 3*8^2 = 481 = numSites
					break;
				case 19:
					return atan<radians>(14.0 / (13.0*sqrt(3.0)));	 // 26^2 + 3*3^2 = 703 = numSites
					break;
			}
			break;
		default:
			error("Invalid cluster size and/or number of tiers");
		}
	}
    error("This point should not be reached");
}

double RegularHexagonalNetwork::compute_tessellation_distance (const size_t numTiers, const size_t clusterSize, const double siteToSiteDistance) const
{
	return sqrt(((numTiers+1)*numTiers*3 + 1.0) * clusterSize) * siteToSiteDistance;
}

void RegularHexagonalNetwork::initialise (const NodeBParams &nodeBParams, const double centreOffset, const radians &startAzimuth)
{
    // create central site and deploy NodeB
    add_clusters(metres(0.0), metres(0.0), nodeBParams, startAzimuth, centreOffset);

    // create surrounding sites and deploy NodeBs at each site
    for (size_t tier = 1; tier <= _numTiers; tier++)
        for (size_t cell = 0; cell < tier; cell++)
        {
            // compute normalised position of site in tier
            const metres x = metres(tier * cos(M_PI/6.0));
	        const metres y = metres(tier * sin(M_PI/6.0) - cell);

	        // convert back to polar co-ordinates
	        metres  r     = sqrt(x*x + y*y);
	        radians theta = atan<radians>(y/x);

	        // translate according to cluster size and cell radius
	        r = r * sqrt(static_cast<double>(_clusterSize)) * siteToSiteDistance;

            switch (_clusterSize)
            {
                case 3:
                case 12:
		            theta -= PI/6.0;
                    break;
		        case 7:
		            theta = theta - atan(sqrt(3.0)/5.0);
                    break;
		        case 13:
		            theta = theta - atan(sqrt(3.0)/7.0);
                    break;
		        case 19:
		            theta = theta - atan(sqrt(3.0)/4.0);
                    break;
                default:
                    break;
	        }

            // rotate for each segment, create site and deploy NodeB
	        for (size_t segment = 1; segment <= 6; segment++)
            {
                theta += PI/3.0;

		        // convert to cartesian co-ordinates
		        const metres x = r * cos(theta);
		        const metres y = r * sin(theta);

                // create site and deploy
                add_clusters(x, y, nodeBParams, startAzimuth, centreOffset);
            }
        }
}

//
// class SquareNetwork
//

SquareNetwork::SquareNetwork (NetworkPropagation &propagation, const metres siteToSiteDistance, const size_t clusterSize, const size_t numTiers, const NodeBParams &nodeBParams, const bool randomiseActiveNodeBs, const double centreOffset, const radians startAzimuth)
    : RegularNetwork (propagation, siteToSiteDistance, M_SQRT1_2 * siteToSiteDistance, clusterSize, numTiers, 9, radians(), ((4*numTiers) + 2) * M_SQRT1_2 * siteToSiteDistance), Network(propagation, randomiseActiveNodeBs)
{
	initialise(nodeBParams, centreOffset, startAzimuth);
}

void SquareNetwork::add_clusters (const metres &x, const metres &y, const NodeBParams &nodeBParams, const radians &startAzimuth)
{
    // only cluster size of unity is valid at present
    Site *site = add_site(x, y, 0, nodeBParams.fadeParams, 0.0);
    add_nodeBs(site, nodeBParams, startAzimuth);
}

void SquareNetwork::initialise (const NodeBParams &nodeBParams, const double centreOffset, const radians &startAzimuth)
{
    const long numTiers = static_cast<long>(_numTiers);

    for (long row = -numTiers; row <= numTiers; row++)
	    for (long column = -numTiers; column <= numTiers; column++)
        {
	        // compute centres of surrounding sites
	        const metres x = siteToSiteDistance * column;
	        const metres y = siteToSiteDistance * row;

            // create site and deploy
            add_clusters(x, y, nodeBParams, startAzimuth);
        }
}

//
// class HexagonalNetwork
//

HexagonalNetwork::HexagonalNetwork (NetworkPropagation &propagation, const metres siteToSiteDistance, const size_t clusterSize, const size_t numTiers, const NodeBParams &nodeBParams, const bool randomiseActiveNodeBs, const double centreOffset, const radians startAzimuth)
	: RegularHexagonalNetwork (propagation, siteToSiteDistance, siteToSiteDistance / sqrt(3.0), clusterSize, numTiers), Network(propagation, randomiseActiveNodeBs)
{
	initialise(nodeBParams, centreOffset, startAzimuth);
}

//
// class CloverLeafNetwork
//

CloverLeafNetwork::CloverLeafNetwork (NetworkPropagation &propagation, const metres siteToSiteDistance, const size_t clusterSize, const size_t numTiers, const NodeBParams &nodeBParams, const bool randomiseActiveNodeBs, const double centreOffset, const radians startAzimuth)
	: RegularHexagonalNetwork (propagation, siteToSiteDistance, 2.0 * siteToSiteDistance / 3.0, clusterSize, numTiers), Network(propagation, randomiseActiveNodeBs)
{
	initialise(nodeBParams, centreOffset, startAzimuth);
}

