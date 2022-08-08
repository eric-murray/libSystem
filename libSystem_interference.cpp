#define _USE_MATH_DEFINES

#include <libEDM_matrix.h>

#include <libSystem_interference.h>
#include <libSystem_networks.h>
#include <libSystem_propagation.h>
#include <libSystem_sites.h>
#include <libSystem_stations.h>

Network *Interference::_network = NULL;

dMatrix Interference::_crossCoupling = dMatrix(0);
dVector Interference::_noise         = dVector(0);

dMatrix Interference::_previousCrossCoupling = dMatrix(0);
dVector Interference::_previousNoise         = dVector(0);

dMatrix Interference::_inverse = dMatrix(0);
dVector Interference::_z       = dVector(0);
dVector Interference::_v       = dVector(0);

dMatrix Interference::compute_crosscoupling_matrix()
{
    // Compute cross-coupling matrix
    dMatrix crossCoupling(_network->numNodeBs(), 0.0);

    for (size_t i1=0; i1<_network->numSites(); i1++)
        for (size_t j1=0; j1<_network->site(i1)->numNodeBs(); j1++)
        {
            NodeB        *nodeB1 = _network->site(i1)->nodeB(j1);
            const size_t  i      = nodeB1->id;

            for (size_t i2=0; i2<_network->numSites(); i2++)
                for (size_t j2=0; j2<_network->site(i2)->numNodeBs(); j2++)
                {
                    NodeB        *nodeB2 = _network->site(i2)->nodeB(j2);
                    const size_t  m      = nodeB2->id;

                    // iterate over UEs connected to NodeB [i1,j1]
                    for (size_t n=0; n<nodeB1->numConnectedPCUEs(); n++)
                    {
                        const PowerControlledUE *ue = nodeB1->connectedPCUEs()->at(n);

						crossCoupling[i][m] += ue->pathgain(m) * ue->crossCouplingMultiplier() * ((i==m) ? orthfactor() : 1.0);
                    }
                }
        }

    return crossCoupling;
}

dVector Interference::compute_noise_vector()
{
    // Compute noise vector
    dVector noise(_network->numNodeBs(), 0.0);
    for (size_t i=0; i<_network->numSites(); i++)
        for (size_t j=0; j<_network->site(i)->numNodeBs(); j++)
            for (size_t n=0; n<_network->site(i)->nodeB(j)->numConnectedPCUEs(); n++)
            {
                const PowerControlledUE *ue      = _network->site(i)->nodeB(j)->powerControlledUE(n);
                const size_t             nodeBId = _network->site(i)->nodeB(j)->id;

                noise[nodeBId] += ue->effectiveNoiseFloor() * ue->crossCouplingMultiplier();
            }

    return noise;
}

bool Interference::compute_nodeb_transmit_powers(const bool check)
{
    // Compute and update NodeB transmit powers
    // Returns true if successful, else false

	const dMatrix _previousInverse = _inverse;

    // Compute inverse of I-K using Sherman-Morrison formula
	const dVector w      = _inverse.transpose() * _v;
	const double  lambda = inner_product(_v, _z);

	_inverse = _inverse - outer_product(_z, w) * (1.0 / (1.0 + lambda));

	if ( check )
	{
		// Test that computed inverse matrix is correct by comparing with directly computed version
		const dMatrix inverse = !(dMatrix(_crossCoupling.size()) - _crossCoupling);
		for (size_t row = 0; row < inverse.size(); row++)
			for (size_t col = 0; col < inverse.size(); col++)
				if ( abs(inverse[row][col] - _inverse[row][col]) > 0.00001 )
					error("Inconsistent inverse matrix");
	}

    // Compute NodeB transmit powers
    const dVector txPowers = _inverse * _noise;

    // Check that computation has been successful
    for (size_t i=0; i<_network->numSites(); i++)
        for (size_t j=0; j<_network->site(i)->numNodeBs(); j++)
        {
            const NodeB  *nodeB        = _network->site(i)->nodeB(j);
            const size_t  nodeBId      = nodeB->id;
            const double  nodeBTxPower = txPowers[nodeBId];
            if ( nodeBTxPower < 0.0 ||
				 nodeBTxPower > nodeB->maxTrafficChannelTxPower ||
				 nodeBTxPower + nodeB->commonChannelTxPower + nodeB->sharedChannelTxPower > nodeB->maxTxPower)
			{
				_inverse       = _previousInverse;
				_crossCoupling = _previousCrossCoupling;
				_noise         = _previousNoise;
                return false;
			}
        }

    // OK, all UEs can be supported, so update NodeB transmit powers
    for (size_t i=0; i<_network->numSites(); i++)
        for (size_t j=0; j<_network->site(i)->numNodeBs(); j++)
        {
                  NodeB  *nodeB   = _network->site(i)->nodeB(j);
            const size_t  nodeBId = nodeB->id;
            nodeB->setTrafficChannelTxPower(txPowers[nodeBId]);
            nodeB->setTxPower(nodeB->commonChannelTxPower + nodeB->sharedChannelTxPower + txPowers[nodeBId]);
        }

	compute_ue_receive_powers();

    return true;
}

void Interference::compute_ue_receive_powers()
{
    // Compute and update UE effective receive powers
    for (size_t i=0; i<_network->numSites(); i++)
		for (size_t j=0; j<_network->site(i)->numNodeBs(); j++)
			for (size_t k=0; k<_network->site(i)->nodeB(j)->connectedUEs()->size(); k++)
			{
				UE *ue = _network->site(i)->nodeB(j)->connectedUEs()->at(k);
				ue->computeEffectiveReceivedPower();
			}
}

void Interference::initialise (Network *network)
{
	_network = network;
}

double Interference::orthfactor()
{
	return _network->propagation.orthfactor;
}

void Interference::recompute()
{
	_crossCoupling = compute_crosscoupling_matrix();
	_noise         = compute_noise_vector        ();
}

void Interference::reset()
{
    _crossCoupling.zeros();
    _noise.zeros();
	_z.zeros();
	_v.zeros();

	// Initialise inverse matrix to identity
	_inverse = dMatrix(_inverse.size());
}

MHz Interference::systemBandwidth()
{
	return _network->propagation.systemBandwidth;
}

bool Interference::update (const PowerControlledUE &ue, const bool check)
{
	if ( _crossCoupling.size() != _network->numNodeBs() )
	{
		// _crossCoupling and _noise have not been correctly initialised, so do so now
		// NOTE - set_size zeroes existing matrix. Need to overload resize so that nodeBs can be added at any time.
	    _crossCoupling.set_size(_network->numNodeBs(), _network->numNodeBs(), 0.0);
        _noise        .resize  (_network->numNodeBs(),                        0.0);

		// recompute inverse matrix
		_inverse = !(dMatrix(_crossCoupling.size()) - _crossCoupling);

		_z.resize(_network->numNodeBs(), 0.0);
		_v.resize(_network->numNodeBs(), 0.0);
	}

    update_crosscoupling_matrix(ue, check);
    update_noise_vector        (ue, check);

	return compute_nodeb_transmit_powers(check);
}

void Interference::update_crosscoupling_matrix (const PowerControlledUE &ue, const bool check)
{
	_previousCrossCoupling = _crossCoupling;
    const size_t i = ue.connectedNodeB()->id;

	// update _z as the negative of the ith column of current inverse
	for (size_t m=0; m<_inverse.size(); m++)
		_z[m] = -_inverse[m][i];

    for (size_t j=0; j<_network->numSites(); j++)
        for (size_t k=0; k<_network->site(j)->numNodeBs(); k++)
        {
            const size_t m = _network->site(j)->nodeB(k)->id;

			// compute mth element of _v
			_v[m] = ue.pathgain(m) * ue.crossCouplingMultiplier() * ((i==m) ? orthfactor() : 1.0);

			// and update cross coupling matrix
			_crossCoupling[i][m] += _v[m];
        }

    if ( check && (_crossCoupling != compute_crosscoupling_matrix()) )
        error("Inconsistent cross-coupling matrix");
}

void Interference::update_noise_vector (const PowerControlledUE &ue, const bool check)
{
	_previousNoise = _noise;
    const size_t i = ue.connectedNodeB()->id;

    _noise[i] += ue.effectiveNoiseFloor() * ue.crossCouplingMultiplier();

    if ( check && (_noise != compute_noise_vector()) )
        error("Inconsistent cross-coupling matrix");
}