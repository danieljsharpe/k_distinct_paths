/*
C++ script to write files "min.data.dummy" and "ts.data.dummy" representing a random kinetic transition network

Compile with:
g++ -std=c++11 random_ktn.cpp -o random_ktn

Usage:
random_ktn <ndim> <xmax> <d_thr> <V> <refp_id> <refp1> <refp2> <errtol> <nbasins> <op> <mu> <sigma> <b> (<seed>)
    <ndim>    - dimensionality of space in which vertices are randomly embedded
    <xmax>    - bounds of space are (-xmax,+xmax) in all dimensions ndim
    <d_thr>   - Euclidean distance threshold for considering two minima to be connected
    <V>       - number of vertices
    <refp_id> - reference degree distribution type. 0 = Poisson, 1 = power law, 2 = Gaussian
    <refp1>   - parameter 1 for ref degree distribution
    <refp2>   - parameter 2 for ref degree distribution (ignored for some distributions)
    <errtol>  - tolerance of error between real and reference degree distributions
    <nbasins> - number of vertices to be randomly selected as harmonic basins
    <op>      - order parameter determining roughness of landscape
    <mu>      - mean of minima energies (Boltzmann factor gives node weights)
    <sigma>   - std dev of minima and transition state energies
    <b>       - TS energy is b+Gaussian(0,sigma) higher in energy than the highest-energy of the
                pair of minima that it connects (also TS energy must be greater)
    <seed>    - (optional) random seed
e.g.
random_ktn 5 5 0.1 1000 0 1 0 0.01 2 1 0 0.5 8

Daniel J. Sharpe
Feb 2019
*/

#include <iostream>
#include <string>
#include <vector>
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

static double rand_unif(double xmin, double xmax, int seed) {
    static std::default_random_engine generator (seed);
    static std::uniform_real_distribution<double> distribution(xmin,xmax);
    return distribution(generator);
}
/*
template <class T>
static T rand_unif(T xmin, T xmax, int seed) {
    static std::default_random_engine generator (seed);
    if (typeid(T) == typeid(double)) {
        static std::uniform_real_distribution<double> distribution(xmin,xmax); }
    else {
        static std::uniform_int_distribution<int> distribution(xmin,xmax); }
    return distribution(generator);
}
rand_unif::distribution(0,1);
*/
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

vector<int> get_deg_distrib(double *edges, const int ndim, const int V) {

    vector<int> deg_distrib(20,0);
    for (int i=0;i<V;i++) {
        int deg = 0;
        for (int j=0;j<V;j++) {
            if (edges[(i*V)+j] != 0.) { deg++; } }
        //cout << "finished with node " << i << " degree is: " << deg << endl;
        try { deg_distrib[deg]++; } catch (const out_of_range& e ) {
            cout << "caught exception, i is: " << i << endl;
            int curr_size = deg_distrib.size();
            for (int k=curr_size+1;k<deg+1;k++) {
                deg_distrib.push_back(0); }
            deg_distrib[deg] = 1;
        }
    }
    return deg_distrib;
}

void get_final_edges (double *v_pos, double *edges, const int ndim, const int V, \
                      double d_thr, double errtol, Degree_distrib_funcs *ddf_obj, \
                      Ddfmembfunc ddf_ptr) {

    double err = numeric_limits<double>::infinity();
    cout << "getting degree distribution..." << endl;
    vector<int> deg_distrib = get_deg_distrib(edges,ndim,V);
    cout << "got degree distribution" << endl;
    int nattempt = 0;
    do {
        err = get_distrib_err(deg_distrib,ddf_obj,ddf_ptr,V);
        if ((nattempt%10000==0)&&(nattempt!=0)) {
            cout << "current error after " << nattempt << " steps: " << err << endl; }
        // Monte Carlo move
        // move a random vertex to a new point in the embedding space

        nattempt++;
    } while (err > errtol and nattempt < 100000);
    if (not nattempt<100000) {
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
                ssqd += pow(v_pos[(i*ndim)+ndim]-v_pos[(j*ndim)+ndim],2.); }
            if (pow(ssqd,0.5) < d_thr) { // edge exists
                edges[(i*V)+j] = 1.; edges[(j*V)+i] = 1.; }
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
            v_pos[(i*j)+j] = rand_unif(-xmax,xmax,seed); }
    }
}

int main(int argc, char** argv) {

    // initialise
    const int ndim = stoi(argv[1]); double xmax = stod(argv[2]);
    double d_thr = stod(argv[3]);
    const int V = stoi(argv[4]); int refp_id = stoi(argv[5]);
    double refp1 = stod(argv[6]); double refp2 = stod(argv[7]);
    double errtol = stod(argv[8]);
    int nbasins = stoi(argv[9]); double op = stod(argv[10]);
    double mu = stod(argv[11]); double sigma = stod(argv[12]);
    double b = stod(argv[13]);
    int seed;
    if (argc > 14) { seed = stoi(argv[14]); }
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
    cout << "found initial edges" << endl;
    get_final_edges(v_pos,edges,ndim,V,d_thr,errtol,degree_distrib_func,ddf_ptr);
    cout << "found final edges" << endl;

    // write networks to files

    // cleanup
    delete degree_distrib_func;
    delete[] v_pos;
    delete[] edges;

    return 0;
}
