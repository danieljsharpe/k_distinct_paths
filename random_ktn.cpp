/*
C++ script to write files "min.data.dummy" and "ts.data.dummy" representing a random kinetic transition network

Compile with:
g++ -std=c++11 random_ktn.cpp -o random_ktn

Usage:
random_ktn <ndim> <xmax> <d_thr> <V> <refp_id> <refp1> <refp2> <maxit> <errtol> <nbasins> <op> <mu> <sigma> <b> (<seed>)
    <ndim>    - dimensionality of space in which vertices are randomly embedded
    <xmax>    - bounds of space are (-xmax,+xmax) in all dimensions ndim
    <d_thr>   - Euclidean distance threshold for considering two minima to be connected
    <V>       - number of vertices
    <refp_id> - reference degree distribution type. 0 = Poisson, 1 = power law, 2 = Gaussian
    <refp1>   - parameter 1 for ref degree distribution
    <refp2>   - parameter 2 for ref degree distribution (ignored for some distributions)
    <maxit>   - max no. of attempts for position changes of vertices
    <errtol>  - tolerance of error between real and reference degree distributions
    <nbasins> - number of vertices to be randomly selected as harmonic basins
                if = -1, then the "pot_func" function is used to assign minima energies, and <op> & <mu>
                are ignored, as is <sigma> for minima energies.
    <op>      - order parameter determining roughness of landscape
    <mu>      - mean of minima energies (Boltzmann factor gives node weights)
    <sigma>   - std dev of minima and transition state energies
    <b>       - TS energy is b+Gaussian(0,sigma) higher in energy than the highest-energy of the
                pair of minima that it connects (also TS energy must be greater)
    <seed>    - (optional) random seed

e.g.1. fit 1000 nodes to a Poisson distribution with lambda=4, with energies according to 2 minima chosen as harmonic basins
./random_ktn 5 5 2.5 1000 0 4 0 1000 0.01 2 1 0 0.5 8 15

e.g.2. fit 1000 nodes to a Poisson distribution with lambda=4, with energies according to the user-defined
       three-hole potential
./random_ktn 2 2 0.15 1000 0 4 0 1000 0.01 -1 0 0 0.2 0.25

random_ktn 5 5 2.5 1000 2 12.4 6.1 10 0.01 2 1 0 0.5 8

Daniel J. Sharpe
Feb 2019
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <random>
#include <typeinfo>
#include <utility>
#include <stdlib.h>
#include <ctime>

using namespace std;

/* a user-defined potential function, called to give minima energies if nbasins = -1 */
double pot_func (double *x) {

    return 3.*exp(-pow(x[0],2.)-pow(x[1]-(1./3.),2.)) - 3.*exp(-pow(x[0],2.)-pow(x[1]-(5./3.),2.)) \
           - 5.*exp(-pow(x[0]-1.,2.)-pow(x[1],2.)) - 5.*exp(-pow(x[0]+1.,2.)-pow(x[1],2.)) \
           + 0.2*pow(x[0],4.) + 0.2*pow(x[1]-(1./3.),4.);
}

struct Degree_distrib_funcs {

    double ddf1, ddf2;
    double pi = acos(-1.0);

    Degree_distrib_funcs(double x, double y) : ddf1(x), ddf2(y) {}

    double poisson_distrib (int k) { // ddf1 = "lambda", ddf2 not used
        return (pow(ddf1,double(k))*exp(-ddf1)/double(tgamma(k+1))); }

    double power_law_distrib (int k) { // alpha*k(-gamma)
        return ddf1*pow(k,-ddf2); }

    double gaussian_distrib (int k) { // ddf1 = mean, ddf2 = std dev
        return (1./(ddf2*pow(2.*pi,0.5)))*exp(-pow(double(k)-ddf1,2.)/(2.*pow(ddf2,2.))); }
};

typedef double (Degree_distrib_funcs::*Ddfmembfunc) (int k);

void write_energies (double *min_ens, double *ts_ens, int V, int n_ts, double *edges, \
                     pair<int,int> *ts_conns, vector<int> *bad_min) {

    bool is_bad;
    ofstream mindataf;
    mindataf.open("min.data.dummy",ios::out);
    for (int i=0;i<V;i++) {
        is_bad = false;
        for (auto &min_id: *bad_min) {
            if (i==min_id) { is_bad = true; } }
        if (is_bad) { continue; }
        mindataf << min_ens[i] << endl; }
    mindataf.close();
    ofstream tsdataf;
    tsdataf.open("ts.data.dummy",ios::out);
    for (int i=0;i<n_ts;i++) {
        tsdataf << ts_ens[i] << "  " << ts_conns[i].first << "  " << ts_conns[i].second << endl; }
   
}

void write_posns (double *v_pos, double *ts_pos, vector<int> *bad_min, int V, int n_ts, int ndim) {

    bool is_bad;
    ofstream minposf;
    minposf.open("min.pos.dummy",ios::out);
    for (int i=0;i<V;i++) {
        is_bad = false;
        for (auto &min_id: *bad_min) {
            if (i==min_id) { is_bad = true; } }
        if (is_bad) { continue; }
        for (int dim=0;dim<ndim;dim++) {
            minposf << v_pos[(i*ndim)+dim] << "   "; }
        minposf << endl;
    }
    minposf.close();
    ofstream tsposf;
    tsposf.open("ts.pos.dummy",ios::out);
    for (int i=0;i<n_ts;i++) {
        for (int dim=0;dim<ndim;dim++) {
            tsposf << ts_pos[(i*ndim)+dim] << "   "; }
        tsposf << endl;
    }
    tsposf.close();
}

double rand_normal(double mean, double std_dev, int seed) {

    static std::default_random_engine generator (seed);
    std::normal_distribution<double> distribution(mean,std_dev);
    return distribution(generator);
}

void get_ts_ens_and_pos (double *min_ens, double *ts_ens, double *ts_pos, int n_ts, double mu, \
                         double sigma, double b, int ndim, double *v_pos, double *edges, int V, \
                         vector<int> *bad_min, pair<int,int> *ts_conns, int seed) {

    int n_bad_i = 0; int n_bad_j; bool is_bad_i; bool is_bad_j; double higher_en;
    int ts_id = 0;
    for (int i=0;i<V;i++) {
        is_bad_i = false;
        for (auto &min_id: *bad_min) {
            if (i==min_id) { is_bad_i = true; } }
        if (is_bad_i) { n_bad_i++; }
        for (int j=0;j<V;j++) {
            n_bad_j = 0; is_bad_j = false;
            for (auto &min_id: *bad_min) {
                if (j==min_id) { is_bad_j = true; } }
            if (is_bad_j) { n_bad_j++; }
            if (edges[(i*V)+j] != 0.) { // we have a transition state
                ts_id++;
                // assign the connections for this transition state
                ts_conns[ts_id] = make_pair(i-n_bad_i,j-n_bad_j);
                // assign an energy to this transition state
                if (min_ens[i] < min_ens[j]) { higher_en = min_ens[i]; }
                else { higher_en = min_ens[j]; }
                do {
                    ts_ens[ts_id] = higher_en + b + rand_normal(mu,sigma,seed);
                    if (ts_ens[ts_id] > higher_en) { break; }
                } while (true);
                // assign a position to this transition state - directly in between the two minima
                for (int dim=0;dim<ndim;dim++) {
                    ts_pos[(ts_id*ndim)+dim] = (v_pos[(i*ndim)+dim] - v_pos[(j*ndim)+dim])/2.; }
            }
        }
    }

}

void get_min_ens (double *min_ens, int V, int nbasins, int op, double mu, double sigma, \
                  double (*pot_func_ptr)(double *), int ndim, double *v_pos) {

    if (pot_func_ptr) { // the user-defined potential function gives the energies of minima
        double init_coords[] = {1.,1.3};
        for (int i=0;i<V;i++) {
            double *pos_i = new double[ndim];
            std::copy(v_pos+(i*ndim),v_pos+(i*ndim)+ndim,pos_i);
            min_ens[i] = pot_func_ptr(pos_i);
            delete[] pos_i;
        }
    } else { // the nbasins and op arguments give the energies of minima
    }
}

template <class T>
static T rand_unif(T xmin, T xmax, int seed) {
    static std::default_random_engine generator (seed);
    if (typeid(T) == typeid(double)) {
        std::uniform_real_distribution<double> distribution1(xmin,xmax);
        return distribution1(generator); }
    else if (typeid(T) == typeid(int)) {
        std::uniform_int_distribution<int> distribution2(xmin,xmax);
        return distribution2(generator); }
}

double get_distrib_err (vector<int> deg_distrib, Degree_distrib_funcs *ddf_obj, \
                        Ddfmembfunc ddf_ptr, const int V) {

    vector<double> deg_distrib_ref(deg_distrib.size(),0.);
    int i =0; double err = 0.; // error between reference and current degree distributions
    for (auto &x: deg_distrib) {
        deg_distrib_ref[i] = (ddf_obj->*ddf_ptr)(i);
        err += (abs(deg_distrib_ref[i]-(double(deg_distrib[i])/double(V))));
        i++;
    }
    return err;
}

vector<int> get_deg_distrib(double *edges, const int ndim, const int V, int old_cap,
                            vector<int> *bad_min=nullptr) {

    vector<int> deg_distrib(old_cap,0);
    vector<int>::size_type sz = deg_distrib.capacity();
    for (int i=0;i<V;i++) {
        int deg = 0;
        for (int j=0;j<V;j++) {
            if (edges[(i*V)+j] != 0.) { deg++; } }
        if (deg > sz) { // update capacity of vector
            deg_distrib.reserve(deg+1);
            for (int k=sz+1;k<deg;k++) {
                deg_distrib.push_back(0);
                sz = deg_distrib.capacity(); }
        }
        if ((bad_min) && (deg==0)) { bad_min->push_back(i); }
        deg_distrib[deg]++;
    }
    return deg_distrib;
}

void update_edges (double *v_pos, double *edges, const int ndim, const int V, \
                   double d_thr) {
}

void get_init_edges (double *, double *, const int, const int, double);

pair<pair<int,int>,vector<int> *> get_final_edges (double *v_pos, double *edges, const int ndim,
                      const int V, \
                      double d_thr, int maxit, double errtol, Degree_distrib_funcs *ddf_obj, \
                      Ddfmembfunc ddf_ptr, double xmax, int seed) {

    int v_id;
    double err = numeric_limits<double>::infinity();
    double old_err; int old_cap = 5; double old_d_thr;
///    vector<int> deg_distrib(old_cap,0);
    vector<int> *bad_min = new vector<int>; // IDs of disconnected minima
    int nattempt = 0;
    do {
        old_err = err;
        if ((nattempt%10==0)&&(nattempt!=0)) { // scale distance threshold move
        old_d_thr = d_thr; double rand_no = rand_unif(0.9,1.05,seed);
        cout << "current error after " << nattempt << " steps: " << err << endl;
        get_init_edges(v_pos,edges,ndim,V,d_thr);
        vector<int> deg_distrib2 = get_deg_distrib(edges,ndim,V,old_cap);
        err = get_distrib_err(deg_distrib2,ddf_obj,ddf_ptr,V);
        if (not (err < old_err)) { // reject move
            cout << "  rejecting threshold move..." << endl;
            err = old_err; old_cap = deg_distrib2.capacity(); d_thr = old_d_thr; }
        else { cout << "  accepting threshold move..." << endl; /// deg_distrib = deg_distrib2;
        }
        } else { // move a random vertex to a new point in the embedding space
        v_id = rand_unif(0,V-1,seed);
        double *old_pos = new double[ndim];
        std::copy(v_pos+(v_id*ndim),v_pos+(v_id*ndim)+ndim,old_pos);
        for (int i=0;i<ndim;i++) {
            v_pos[(v_id*ndim)+i] = rand_unif(-xmax,xmax,seed); }
        // update degree distribution and accept if better match with reference distribution
        get_init_edges(v_pos,edges,ndim,V,d_thr);
////        update_edges(v_pos,edges,ndim,V,d_thr);
//        cout << "getting deg_distrib2..." << endl;
        vector<int> deg_distrib2 = get_deg_distrib(edges,ndim,V,old_cap);
//        cout << "got deg_distrib2..." << endl;
        err = get_distrib_err(deg_distrib2,ddf_obj,ddf_ptr,V);
//        cout << "got distribution error..." << endl;
        old_cap = deg_distrib2.capacity();
        if (not (err < old_err)) { // reject move
//            cout << "rejecting move..." << endl;
            err = old_err;
            for (int i=0;i<ndim;i++) { v_pos[(v_id*ndim)+i] = old_pos[i]; }
        } else { cout << "accepting move..." << endl;
///                 deg_distrib.reserve(deg_distrib2.capacity());
///                 deg_distrib = deg_distrib2;
        }
//        cout << "deleting old_pos..." << endl;
        delete[] old_pos;
//        cout << "deleted old_pos..." << endl;
        }
        nattempt++;
        if ((err < errtol) || (nattempt > maxit)) { break; }
    } while (true);
    if (not nattempt<maxit) {
        cout << "Warning: exceeded max number of attempts to adjust degree distribution" << endl;
        cout << "error between reference and current distributions: " << err << endl; }
    else { cout << "succeeded to decrease error below tolerance, final error: " << err << endl; }
    cout << "final degree distribution (vs ref):" << endl;
    vector<int> deg_distrib(old_cap,0);
    deg_distrib = get_deg_distrib(edges,ndim,V,old_cap,bad_min);
    int i=0; int n_ts=0;
    for (auto &x: deg_distrib) {
        n_ts += i*deg_distrib[i];
        cout << "  " << i << "   " << deg_distrib[i] << "      " << \
                double(deg_distrib[i])/double(V) << "   " << (ddf_obj->*ddf_ptr)(i) << endl; i++; }
    pair<int,int> retval_tmp = make_pair(deg_distrib[0],n_ts);
    pair<pair<int,int>,vector<int> *> retval = make_pair(retval_tmp,bad_min);
    return retval;
}

void get_init_edges (double *v_pos, double *edges, const int ndim, const int V, \
                     double d_thr) {

    for (int i=0;i<V;i++) {
        for (int j=i+1;j<V;j++) {
            double ssqd = 0.; // sum of square distances
            for (int dim=0;dim<ndim;dim++) {
                ssqd += pow(v_pos[(i*ndim)+dim]-v_pos[(j*ndim)+dim],2.); }
            if (pow(ssqd,0.5) < d_thr) { // edge exists
                edges[(i*V)+j] = 1.; edges[(j*V)+i] = 1.; }
            else { edges[(i*V)+j] = 0.; edges[(j*V)+i] = 0.; } // edge does not exist
        }
    }
}

/*
// works if dimensions of arrays to be passed known at compile time
template <size_t V, size_t ndim>
void get_v_pos (double (&v_pos)[V][ndim]) {}
*/

void get_v_pos (double *v_pos, const int ndim, const int V, double xmax, int seed) {

    for (int i=0;i<V;i++) {
        for (int j=0;j<ndim;j++) {
            v_pos[(i*ndim)+j] = rand_unif(-xmax,xmax,seed); }
    }
}

int main(int argc, char** argv) {

    // initialise
    const int ndim = stoi(argv[1]); double xmax = stod(argv[2]);
    double d_thr = stod(argv[3]);
    const int V = stoi(argv[4]); int refp_id = stoi(argv[5]);
    double refp1 = stod(argv[6]); double refp2 = stod(argv[7]);
    int maxit = stoi(argv[8]); double errtol = stod(argv[9]);
    int nbasins = stoi(argv[10]); double op = stod(argv[11]);
    double mu = stod(argv[12]); double sigma = stod(argv[13]);
    double b = stod(argv[14]);
    int seed;
    if (argc > 15) { seed = stoi(argv[15]); }
    else { srand(time(NULL)); seed = rand(); }
    Degree_distrib_funcs *degree_distrib_func = new Degree_distrib_funcs(refp1,refp2);
    Ddfmembfunc ddf_ptr;
    if (refp_id==0) { ddf_ptr = &Degree_distrib_funcs::poisson_distrib; }
    else if (refp_id==1) { ddf_ptr = &Degree_distrib_funcs::power_law_distrib; }
    else if (refp_id==2) { ddf_ptr = &Degree_distrib_funcs::gaussian_distrib; }
    double (*pot_func_ptr)(double *);
    if (nbasins==-1) { pot_func_ptr = &pot_func; }
    else { pot_func_ptr = nullptr; }
    cout << "finished initialising" << endl;

    // build kinetic transition network connectivity from predefined degree distribution
    double *v_pos = new double[V*ndim]; // array of vertex positions
    double *edges = new double[V*V]; // array of edge energies (=0. if edge nonexistent)
    get_v_pos(v_pos,ndim,V,xmax,seed);
    get_init_edges(v_pos,edges,ndim,V,d_thr);
    cout << "assigned initial edges" << endl;
    pair<pair<int,int>,vector<int> *> gfe_ret = get_final_edges(v_pos,edges,ndim,V,d_thr,maxit, \
                    errtol,degree_distrib_func,ddf_ptr,xmax,seed);
    pair<int,int> n_sp = gfe_ret.first;
    vector<int> *bad_min = gfe_ret.second;
    cout << "assigned final edges" << endl;
    cout << "there are " << n_sp.first << " minima that are disconnected and will be ignored" << endl;

    // assign minima and transition state energies
    int n_minima = V - n_sp.first; // number of connected minima
    int n_ts = n_sp.second; // number of transition states
    cout << "there are " << n_ts << " transition states" << endl;
    double *min_ens = new double[V];
    double *ts_ens = new double[n_ts];
    double *ts_pos = new double[n_ts*ndim];
    pair<int,int> *ts_conns = new pair<int,int>[n_ts];
    get_min_ens(min_ens,V,nbasins,op,mu,sigma,pot_func_ptr,ndim,v_pos);
    get_ts_ens_and_pos(min_ens,ts_ens,ts_pos,n_ts,mu,sigma,b,ndim,v_pos,edges,V,bad_min, \
                       ts_conns,seed);
    cout << "assigned minima and transition state energies" << endl;

    // write kinetic transition network to files
    write_energies(min_ens,ts_ens,V,n_ts,edges,ts_conns,bad_min);
    write_posns(v_pos,ts_pos,bad_min,V,n_ts,ndim);
    cout << "finished writing min.data and ts.data files" << endl;

    // cleanup
    delete degree_distrib_func;
    delete[] v_pos;
    delete[] edges;
    delete[] min_ens; delete[] ts_ens; delete[] ts_pos; delete[] ts_conns;
    delete[] bad_min;

    return 0;
}
