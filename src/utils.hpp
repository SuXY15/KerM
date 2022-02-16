#include <math.h>
#include <vector>
#include <string>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>

#ifdef DOUBLE
    typedef double real;
#else
    typedef float real;
#endif
typedef std::vector<int> ivec;
typedef std::vector<real> vec;
typedef std::chrono::high_resolution_clock timer;
typedef std::chrono::duration<float, std::milli> duration;

#define PI 3.141592653
#define REAL_MIN -1e10
#define REAL_MAX  1e10

#define min(a,b) (a)>(b)?(b):(a)
#define max(a,b) (a)>(b)?(a):(b)

using namespace std;

class pdfFunc
{
public:
    real xmin, xmax;
public:
    pdfFunc(){};
    pdfFunc(real* xlim):xmin(xlim[0]), xmax(xlim[1]){};
    ~pdfFunc(){};
    real eval(real x) const;
};

void raiseNotDefined()
{
    cerr << "Raise Not Defined." << endl;
}

// srand( (unsigned)time( NULL ) );
real randr()
{
    return (double) rand() / (RAND_MAX);
}

inline void get_minmax(const vector<real> &W, real &Wmin, real &Wmax)
{
    Wmin = REAL_MAX, Wmax = REAL_MIN;
    for(int i=0; i<W.size(); ++i){
        Wmin = min(W[i], Wmin);
        Wmax = max(W[i], Wmax);
    }
}

inline real mean(const vector<real> &data)
{
    real Dsum = 0;
    for(int i=0; i<data.size(); ++i){
        Dsum += data[i];
    }
    return Dsum / data.size();
}

inline real mean(const vector<real> &data, const vector<real> &weights)
{
    real Wsum = 0, Dsum = 0;
    for(int i=0; i<data.size(); ++i){
        Wsum += weights[i];
        Dsum += weights[i] * data[i];
    }
    return Dsum / Wsum;
}

inline real meanAB(const vector<real> &A, const vector<real> &B, const vector<real> &weights)
{
    real Wsum = 0, ABsum = 0;
    for(int i=0; i<weights.size(); ++i){
        Wsum += weights[i];
        ABsum += weights[i] * A[i] * B[i];
    }
    return ABsum / Wsum;
}

inline real variance(const vector<real> &data, const vector<real> &weights)
{
    real Wsum = 0, Dvar = 0;
    real Davr = mean(data, weights);
    for(int i=0; i<data.size(); ++i){
        Wsum += weights[i];
        Dvar += weights[i] * pow(data[i] - Davr, 2);
    }
    return Dvar / Wsum;
}

void readXY(string filename, vector<real> &X, vector<real> &Y)
{
    X.resize(0);
    Y.resize(0);
    real x, y;
    ifstream f(filename, ios::in);
    while(!f.eof()){
        f >> x >> y;
        X.push_back(x);
        Y.push_back(y);
    }
    f.close();
    return;
}

void saveXY(string filename, const vector<real> X, const vector<real> Y)
{
    ofstream f(filename, ios::out);
    f << std::scientific << setprecision(8);
    for(int i=0; i<X.size(); ++i){
        f << setw(15) << X[i] << " " 
          << setw(15) << Y[i] << endl;
    }
    f.close();
    return;
}

void saveCost(string filename, const int N, const real timecost)
{
    ofstream f(filename, ios_base::app);
    f << std::scientific << setprecision(8);
    f << setw(15) << N << " " << setw(15) << timecost << endl;
    f.close();
    return;
}

//---------------------------------------------------------
// FUNC: use binary search to find a location
template <typename T>
inline int binarySearch(int l, int r, T val, vector<T> arr)
{
    int m;
    while (r>l) {
        m = (l + r) >> 1;
        arr[m] <= val ? r = m : l = m + 1;
    }
    return l;
}

//---------------------------------------------------------
// FUNC: quick sort with assistance array
template <typename T1, typename T2>
void quickSort(int l, int r, vector<T1> &A, vector<T2> &nA)
{
    int i, j;
    T1 t1, t2;
    T2 nt1, nt2;
    if (l>r) return;
    t2 = A[l];
    i = l; j = r;
    while (i != j) {
        while (A[j] <= t2 && i<j) j--;
        while (A[i] >= t2 && i<j) i++;
        if (i<j) {
            t1 = A[i];    A[i] = A[j];   A[j] = t1;
            nt1 = nA[i]; nA[i] = nA[j]; nA[j] = nt1;
        }
    }
    A[l] = A[i];  A[i] = t2;
    nt2 = nA[l]; nA[l] = nA[i]; nA[i] = nt2;
    quickSort(l, i - 1, A, nA);
    quickSort(i + 1, r, A, nA);
}

inline void setWeights(vector<real> &particles, vector<real> &weights,
                       vector<real> &prob, string type="uniform")
{
    int N = particles.size(); // N: sample numbers
    weights.resize(N);
    vector<int> index(N);

    for(int i=0; i<N; ++i){
        index[i] = i;
    }

    if (type=="uniform"){
        for(int i=0; i<N; ++i){
            weights[i] = 1./N;
        }
    }
    else if(type=="weighted"){
        int j, l, r;
        real width, Wsum = 0;
        vector<real> v = particles;
        quickSort(0, N-1, v, index);
        for(int i=0; i<N; ++i){
            l = max(i-1, 0);
            r = min(i+1, N-1);
            j = index[i];
            width = (v[r] - v[l])/(r-l);
            weights[j] = prob[j]*width;
            Wsum += weights[j];
        }
        for(int i=0; i<N; ++i){
            weights[i] /= Wsum;
        }

        // cout << "Generating Weighted Samples" << endl;
        // for(int i=0; i<N; ++i){
        //     printf("%7.3e %7.3e %7.3e\n", particles[i], prob[i], weights[i]);
        // }
    }else{
        cout << "Only 'uniform' and 'weighted' are supported." << endl;
        return;
    }
    return;
}

inline void acceptSampling(const pdfFunc f, real* xlim, int N,
                           vector<real> &particles, vector<real> &weights,
                           string type="uniform")
{
    real x, p;
    particles.resize(N);
    vector<real> prob(N);

    for(int i=0; i<N; ++i){
        while(true){
            x = randr()*(xlim[1]-xlim[0]) + xlim[0];
            p = f.eval(x);
            if(randr() <= p) break;
        }
        particles[i] = x;
        prob[i] = p;
    }
    setWeights(particles, weights, prob, type);
    return;
}

inline void acceptSampling(vector<real> &phi, vector<real> &P, int N,
                           vector<real> &particles, vector<real> &weights,
                           string type="uniform")
{
    int k, kp, n = phi.size();
    quickSort(0, n-1, phi, P);
    real Pmin, Pmax;
    get_minmax(P, Pmin, Pmax);

    real x, p;
    real xlim[2] = {phi[0], phi[n-1]};
    particles.resize(N);
    vector<real> prob(N);

    for(int i=0; i<N; ++i){
        while(true){
            x = randr()*(xlim[1]-xlim[0]) + xlim[0];
            k = binarySearch(0, n-1, x, phi); // k <= n-2
            kp = k+1;
            while(kp<n-1 && phi[k+1]-phi[k]<=0) kp++;
            p = ((P[kp]-P[k])/(phi[kp]-phi[k])*(x-phi[k]) + P[k])/Pmax;
            if(randr() <= p) break;
        }
        // printf("%3d %8.3e %8.3e %8.3e %8.3e %8.3e\n",
        //         k, x, p, phi[k], phi[k+1], P[k]);
        particles[i] = x;
        prob[i] = p;
    }
    setWeights(particles, weights, prob, type);
    return;
}
