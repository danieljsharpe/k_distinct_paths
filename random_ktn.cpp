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
    <op>      - order parameter determining roughness of landscape
    <mu>      - mean of minima energies (Boltzmann factor gives node weights)
    <sigma>   - std dev of minima and transition state energies
    <b>       - TS energy is b+Gaussian(0,sigma) higher in energy than the highest-energy of the
                pair of minima that it connects (also TS energy must be greater)
    <seed>    - (optional) random seed
e.g.
random_ktn 5 5 2.5 1000 0 1 0 10 0.01 2 1 0 0.5 8

Daniel J. Sharpe
Feb 2019
*/

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <random>
#include <typeinfo>
#include <stdlib.h>
#include <ctime>

using namespace std;

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

template <class T>
static T rand_unif(T xmin, T xmax, int seed) {
    static std::default_random_engine generator (seed);
    if (typeid(T) == typeid(double)) {
        static std::uniform_real_distribution<double> distribution1(xmin,xmax);
        return distribution1(generator); }
    else if (typeid(T) == typeid(int)) {
        static std::uniform_int_distribution<int> distribution2(xmin,xmax);
        return distribution2(generator); }
}

double get_distrib_err (vector<int> deg_distrib, Degree_distrib_funcs *ddf_obj, \
                        Ddfmembfunc ddf_ptr, const int V) {

    vector<double> deg_distrib_ref(deg_distrib.size(),0.);
    int i =0; double err = 0.; // % error between reference and current degree distributions
    //cout << "getting distribn error..." << endl;
    for (auto &x: deg_distrib) {
        //cout << "k: " << i << " number of nodes: " << deg_distrib[i] << endl;
        deg_distrib_ref[i] = (ddf_obj->*ddf_ptr)(i);
        err += (abs(deg_distrib_ref[i]-(double(deg_distrib[i])/double(V)))/double(V))*100.;
        i++;
    }
    return err;
}

vector<int> get_deg_distrib(double *edges, const int ndim, const int V, int old_cap) {

    vector<int> deg_distrib(old_cap,0);
    vector<int>::size_type sz = deg_distrib.capacity();
    cout << "beginning loop..." << endl;
    for (int i=0;i<V;i++) {
        int deg = 0;
        for (int j=0;j<V;j++) {
            if (edges[(i*V)+j] != 0.) { deg++; } }
        //cout << "finished with node " << i << " degree is: " << deg << endl;
        if (deg > sz) { // update capacity of vector
            cout << "updating capacity of vector... deg = " << deg << endl;
            deg_distrib.reserve(deg+1);
            for (int k=sz+1;k<deg;k++) {
                deg_distrib.push_back(0);
                sz = deg_distrib.capacity();
                cout << "updated capacity of vector to: " << sz << " because of node of deg: " << deg << endl; }
        }
        deg_distrib[deg]++;
    }
    cout << "returning deg_distrib" << endl;
    return deg_distrib;
}

void update_edges (double *v_pos, double *edges, const int ndim, const int V, \
                   double d_thr) {
}

void get_init_edges (double *, double *, const int, const int, double);

void get_final_edges (double *v_pos, double *edges, const int ndim, const int V, \
                      double d_thr, int maxit, double errtol, Degree_distrib_funcs *ddf_obj, \
                      Ddfmembfunc ddf_ptr, double xmax, int seed) {

    double err = numeric_limits<double>::infinity();
    double old_err;
    vector<int> deg_distrib; deg_distrib.reserve(20);
    int nattempt = 0;
    do {
        if ((nattempt%10==0)&&(nattempt!=0)) {
            cout << "current error after " << nattempt << " steps: " << err << endl; }
        old_err = err;
        // Monte Carlo move
        // move a random vertex to a new point in the embedding space
        int v_id = rand_unif(0,V,seed);
        cout << "selected node: " << v_id << endl;
        double *old_pos = new double[ndim];
        std::copy(v_pos+(v_id*ndim),v_pos+(v_id*ndim)+ndim,old_pos);
        for (int i=0;i<ndim;i++) {
            v_pos[(v_id*ndim)+i] = rand_unif(-xmax,xmax,seed); }
        // update degree distribution and accept if better match with reference distribution
        get_init_edges(v_pos,edges,ndim,V,d_thr);
        cout << "assigning deg_distrib..." << endl;
        deg_distrib = get_deg_distrib(edges,ndim,V,deg_distrib.capacity());
        cout << "finished assigning deg_distrib..." << endl;
        err = get_distrib_err(deg_distrib,ddf_obj,ddf_ptr,V);
        cout << "new error: " << err << endl;
        if (not (err < old_err)) { // reject move
            err = old_err;
            std::copy(v_pos+(v_id*ndim),v_pos+(v_id*ndim)+ndim,old_pos); } // restore old position
        delete[] old_pos;
        nattempt++;
    } while (err > errtol and nattempt < maxit);
    if (not nattempt<maxit) {
        cout << "Warning: exceeded max number of attempts to adjust degree distribution" << endl;
        cout << "error between reference and current distributions: " << err << endl; }
    else { cout << "succeeded to decrease error below tolerance, final error: " << err << endl; }
    cout << "final degree distribution:" << endl;
    int i=0;
    for (auto &x: deg_distrib) {
        cout << "  " << i << "   " << deg_distrib[i] << endl; i++; }
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
    cout << "finished initialising" << endl;

    // build kinetic transition network
    double *v_pos = new double[V*ndim]; // array of vertex positions
    double *edges = new double[V*V]; // array of edge energies (=0. if edge nonexistent)
    get_v_pos(v_pos,ndim,V,xmax,seed);
    get_init_edges(v_pos,edges,ndim,V,d_thr);
    cout << "assigned initial edges" << endl;
    get_final_edges(v_pos,edges,ndim,V,d_thr,maxit,errtol,degree_distrib_func, \
                    ddf_ptr,xmax,seed);
    cout << "assigned final edges" << endl;

    cout << "assigned minima and transition state energies" << endl;

    // write kinetic transition network to files

    // cleanup
    delete degree_distrib_func;
    delete[] v_pos;
    delete[] edges;

    return 0;
}
