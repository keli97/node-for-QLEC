#include "simulator.h"

int main(int argc, const char * argv[])
{
    (new Simulator);
    Simulator &sim = Simulator::instance();
    sim.init();
    sim.run();

    //cout<<100<<endl;

    /*for(int i = 0; i<NODE_NUMBER; i++)
    {
        for(int j = 0; j<NODE_NUMBER; j++)
        {
            if(j == NODE_NUMBER-1)
                cout<<sim.P_ij[i][j]<<endl;
            else cout<<sim.P_ij[i][j]<< " ";
        }

    }*/
    if (sim.totalPacketNum > 0)
    {
        cout << "Packet delivery ratio: " << sim.acceptPacketNum / double(sim.totalPacketNum) << endl;
        cout << "Latency per Packet: "<< sim.totalLatency /sim.acceptPacketNum<<endl;
        //cout<< "Energy Cost: "<<sim.consume_energy<<endl;
    }
    else
    {
        cout << "Total packet number is 0!!!" << endl;
    }

    for (int i = 0; i < NODE_NUMBER-1; ++i)
    {
        cout << "Rest energy of node " << i << ": " << sim.nodes_[i]->energy() << endl;
    }

    return 0;
}
