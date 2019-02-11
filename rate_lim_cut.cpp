/*
C++ script to find the rate-limiting cut of a kinetic transition network

Compile with:
g++ -std=c++11 rate_lim_cut.cpp -o rate_lim_cut

Usage:
rate_lim_cut <nmin> <nts> <ts.data file>

Output:
A file "rate_lim_cut.dat" containing the transition states constituting the rate-limiting cut,
in the format: ts_energy / ts_id, listed in order of increasing energy.

See:
Sharpe, D. J. and Wales, D. J.... (in preparation)
Noe, F. et al., J. Chem. Theory Comput., 2006, 2. 840-857

Daniel J. Sharpe
Feb 2019
*/

#include <iostream>
#include <string>
#include <utility>
//#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/stoer_wagner_min_cut.hpp>

using namespace std;

typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::bidirectionalS, 
            boost::no_property,boost::property<boost::edge_weight_t,double>> UndirectedGraph;
typedef boost::property_map<UndirectedGraph,boost::edge_weight_t>::type weight_map_type;
typedef boost::property_traits<weight_map_type>::value_type weight_type;
//typedef pair<double *,int **> ts_info;
typedef pair<double *, int (*)[2]> ts_info;

struct edge_type {
    unsigned int first;
    unsigned int second;
};

ts_info read_ts_data(const int nmin, const int nts, string tsdataf) {

    double lowest_ts_en;
    //double ts_wts[nts] = {1.,2.,3.};
    //double *ts_wts_ptr = &ts_wts[0];
    double *ts_wts = new double[nts];
    ts_wts[0] = 1.; ts_wts[1] = 2.; ts_wts[2] = 3.;
    /*
    int **ts_conns = new int*[nts];
    for (int i=0;i<nts;i++) {
        ts_conns[i] = new int[2];
        ts_conns[i][0] = 5; ts_conns[i][1] = 6;
        cout << i << "  " << ts_conns[i][0] << "  " << ts_conns[i][1] << endl;
    }
    */
    int (*ts_conns)[2] = new int[nts][2];
    for (int i=0;i<nts;i++) {
        ts_conns[i][0] = 5; ts_conns[i][1] = 6; }

    return make_pair(ts_wts, ts_conns);
}


int main(int argc, char *argv[]) {

    const int nmin = stoi(argv[1]);
    const int nts = stoi(argv[2]);
    string tsdataf = argv[3];

    ts_info ts_info1 = read_ts_data(nmin,nts,tsdataf);

    double *ts_wts_ptr = ts_info1.first;
    //int **ts_conns_ptr = ts_info1.second;
    int (*ts_conns_ptr)[2] = ts_info1.second;

    for (int i=0;i<nts;i++) {
        cout << i << "  " << *(ts_wts_ptr+i) << endl;
    }
    for (int i=0;i<nts;i++) {
        cout << i << "  " << ts_conns_ptr[i][0] << "  " << ts_conns_ptr[i][1] << endl;
    }

    // construct the boost graph representing the kinetic transition network
    weight_type *ktn_wts = ts_wts_ptr;
    edge_type ktn_edges[nts];
    for (int i=0;i<nts;i++) {
        ktn_edges[i].first = ts_conns_ptr[i][0]; ktn_edges[i].second = ts_conns_ptr[i][1]; }
    UndirectedGraph ktn(ktn_edges,ktn_edges+nts,ktn_wts,nmin,nts);

    // run the Stoer-Wagner algorithm from the boost library
    //BOOST_AUTO ...

    // cleanup
    delete[] ts_wts_ptr;
    /*
    for (int i=0;i<nts;i++) {
        delete[] ts_conns_ptr[i]; }
    delete [] ts_conns_ptr;
    */
    delete[] ts_conns_ptr;

    return 0;
}
