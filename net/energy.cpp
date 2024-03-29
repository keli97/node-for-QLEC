#include "energy.h"
#include "simulator.h"


const double Energy :: energyPerBitTx = ENERGY_PERBITTX;
const double Energy :: energyPerBitRx = ENERGY_PERBITRX;

// Reduce energy level according operation type. Then check whether it is required to broadcast a relaxation packet.
void Energy :: spend(Node *n, Node*m) {
    Simulator &sim = Simulator::instance();
    if(n == sim.bs_)
        return;
    double energy = n->energy();
    double consume = 0;
    double dis_k;
    dis_k = sim.dis_[n->id()][m->id()];
    double d_thres;
    d_thres = sqrt(EPISILON_FS/EPISILON_MP);
    if(n->nodeType() == TCH)
    {
        if(dis_k<d_thres)
        {
            consume = PACKET_SIZE*COMPRESS_RATIO*EPISILON_FS*pow(dis_k, 2);
        }
        else consume = PACKET_SIZE*COMPRESS_RATIO*EPISILON_MP*pow(dis_k,4);
    }
    else
    {
        if(dis_k<d_thres)
        {
            consume = PACKET_SIZE*EPISILON_FS*pow(dis_k, 2);
        }
        else consume = PACKET_SIZE*EPISILON_MP*pow(dis_k,4);
    }

    pthread_mutex_lock(&sim.mutex_energy);
    sim.energy += consume;
    pthread_mutex_unlock(&sim.mutex_energy);

    n->energy(energy - consume);
    cout<<n->id()<<" consume: "<< consume<<endl;
    cout<<n->id()<<" to "<<m->id()<<endl;
    // Check for threshold
    if(n->reachedThreshold() == -1 && sim.is_alive[n->id()]==true)
    {
        sim.killNode(n);
    }
    if(n->energy()<0)
        n->energy(0);
}
