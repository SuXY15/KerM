#include "utils.hpp"

// base class
class MixingModel
{
public:
    int N;
    real var;
    string name;
    vector<real> phis;
    vector<real> weights;

    timer::time_point t0;
    timer::time_point t1;
    real timecost;
    
public:
    MixingModel(){}

    MixingModel(const vector<real> &particles, const vector<real> &weights):
        name("MixingModel"),
        N(particles.size()),
        phis(particles),
        weights(weights),
        var(variance(particles, weights)),
        timecost(0.0)
    {}

    ~MixingModel(){}
    
    void tic(){
        t0 = timer::now();
    }

    void toc(){
        t1 = timer::now();
        timecost += (duration(t1-t0).count())/1e3;
    }

    virtual void update(real Omega_phi, real dt){
        raiseNotDefined();
    }
};

// IEM model
class IEM: public MixingModel
{
public:
    IEM(){}

    IEM(const vector<real> &particles, const vector<real> &weights):
        MixingModel(particles, weights)
    {
        this->name = "IEM";
    }

    void update(real Omega_phi, real dt)
    {
        tic();
        real omdt = Omega_phi * dt;
        real phi_avr = mean(phis, weights);
        for(int i=0; i<N; ++i){
            phis[i] += -1./2. * omdt * (phis[i] - phi_avr);
        }
        toc();
    }
};

// Modified Curl Model
class MC: public MixingModel
{
public:
    MC(){}

    MC(const vector<real> &particles, const vector<real> &weights):
        MixingModel(particles, weights)
    {
        this->name = "MC";
    }

    void update(real Omega_phi, real dt)
    {
        tic();
        int i, j, p, q, nmix;
        real num_mix, pmix, Wmin, Wmax, Wavr, alpha, phi_pq;
        
        // set number of mixing pairs
        num_mix = 3./2. * N * Omega_phi * dt;
        nmix = floor( num_mix + 1);
        pmix = num_mix / nmix;

        Wavr = mean(weights);
        get_minmax(weights, Wmin, Wmax);

        // mixing
        real* W = weights.data();
        for(int imix=0; imix<nmix; ++imix){
            // select p: marginal prob ~ Wavr + W[p]
            do {
                p = floor(randr()*N);
            } while (W[p] + Wavr <= randr() * (Wmax + Wavr));
            // select q: conditional prob ~ W[p] + W[q]
            do {
                q = floor(randr()*N);
            } while (W[q] + W[p] <= randr() * (Wmax + W[p]));
            if (p==q) continue;

            if (randr()<=pmix){
                alpha = randr();
                phi_pq = (phis[p]*W[p] + phis[q]*W[q])/(W[p]+W[q]);
                phis[p] = phis[p] - alpha * (phis[p] - phi_pq);
                phis[q] = phis[q] - alpha * (phis[q] - phi_pq);
            }
        }
        toc();
    }
};

// EMST Model
class EMST: public MixingModel
{
public:
    vector<real> Bv, dphi;
    vector<int> index;

public:
    EMST(){}

    EMST(const vector<real> &particles, const vector<real> &weights):
        MixingModel(particles, weights)
    {
        this->name = "EMST-1D";
        Bv.resize(N);
        dphi.resize(N);
        index.resize(N);
    }

    void update(real Omega_phi, real dt)
    {
        tic();
        int mv, nv;
        real Wi, wv;
        real A, B, C, alpha, dt_in = dt;

        // inner loop for smaller dt
        while(dt_in>0){
            vector<real> v = phis;

            for(int i=0; i<N; ++i){
                dphi[i] = 0.;
                index[i] = i;
            }
            quickSort(0, N-1, v, index);
            
            Wi = 0;
            for(int i=0; i<N; ++i){
                Wi += weights[index[i]];
                wv = min(Wi, 1-Wi);
                Bv[i] = 2*wv;
            }
            for(int v=0; v<N-1; ++v){
                mv = index[v];
                nv = index[v+1];
                dphi[mv] += - Bv[v]*(phis[mv]-phis[nv])/weights[mv];
                dphi[nv] += - Bv[v]*(phis[nv]-phis[mv])/weights[nv];
            }
            A = meanAB(dphi, dphi, weights);
            B = meanAB(dphi, phis, weights) * 2;
            C = Omega_phi * variance(phis, weights);

            dt = 1.0 * B*B / (4*A*C);
            A = A * dt * dt;
            B = B * dt;
            C = C * dt;
            alpha = (-B+sqrt(abs(B*B-4*A*C)))/2/A;

            dt = min(dt, dt_in);
            // cout << "EMST dt_in=" << dt_in << ", dt=" << dt << endl;

            for(int i=0; i<N; ++i){
                phis[i] += dphi[i] * alpha * dt;
            }
            dt_in -= dt;
        }
        toc();
    }
};

// Kernel Mixing Model
class KerM: public MixingModel
{
public:
    int Nbin;
    real sigma_k;
    real sigma_k2m4;
    vector<int> index;
    vector<real> cdf;
    vector<int> binCount;
    vector<int> binAccum;
    vector<real> binVmin;
    vector<real> binVmax;

public:
    KerM(){}

    KerM(const vector<real> &particles, const vector<real> &weights,
         const real sigma_k = 100., const int Nbin=50):
         MixingModel(particles, weights)
    {
        this->name = "KerM";
        this->sigma_k = sigma_k;
        this->sigma_k2m4 = 4*sigma_k*sigma_k;
        this->cdf.resize(N);
        this->index.resize(N);
        this->Nbin = min(max(Nbin,10), int(sqrt(N)));
        this->binCount.resize(Nbin);
        this->binAccum.resize(Nbin);
        this->binVmax.resize(Nbin);
        this->binVmin.resize(Nbin);
    }

    double kernel_func(real d)
    {
        return exp(-d*d/sigma_k2m4);
    }

    void update(real Omega_phi, real dt)
    {
        tic();
        int i, j, k, p, q, nmix;
        real coeff, f_pq, dvar = 0;
        real Wavr, Wmin, Wmax, Wi = 0;

        // // quicksort for CDF
        // vector<real> v = phis;
        // for(int i=0; i<N; ++i){
        //     index[i] = i;
        // }
        // quickSort(0, N-1, v, index);
        // for(int i=0; i<N; ++i){
        //     Wi += weights[index[i]];
        //     cdf[index[i]] = Wi;
        // }

        // bucketsort for CDF
        real vmin, vmax, vavr;
        get_minmax(phis, vmin, vmax);
        vmin = vmin - 1e-8;
        vmax = vmax + 1e-8;
        vavr = (vmax-vmin) / Nbin;
        for(int k=0; k<Nbin; ++k){
            binCount[k] = 0;
            binVmin[k] = vavr*(k  ) + vmin;
            binVmax[k] = vavr*(k+1) + vmin;
        }
        for(int i=0; i<N; ++i){
            k = floor((phis[i] - vmin) / vavr);
            binCount[k] ++;
        }
        int accum = 0;
        for(int k=0; k<Nbin; ++k){
            binAccum[k] = accum;
            accum += binCount[k];
        }
        for(int i=0; i<N; ++i){
            k = floor((phis[i] - vmin) / vavr);
            cdf[i] = ((phis[i]-binVmin[k])/vavr*binCount[k]+binAccum[k])/N;
        }
        
        Wavr = mean(weights);
        get_minmax(weights, Wmin, Wmax);

        // get coeffs by random selecting N pairs
        real* W = weights.data();
        for(int i=0; i<N; ++i){
            // select p: marginal prob ~ Wavr + W[p]
            do {
                p = floor(randr()*N);
            } while (W[p] + Wavr <= randr() * (Wmax + Wavr));
            // select q: conditional prob ~ W[p] + W[q]
            do {
                q = floor(randr()*N);
            } while (W[q] + W[p] <= randr() * (Wmax + W[p]));
            f_pq = kernel_func(cdf[p]-cdf[q]);
            dvar += f_pq * W[p]*W[q]/(W[p]+W[q])*pow(phis[p]-phis[q],2);
        }
        coeff = variance(phis, weights) / dvar;

        // // get coeffs by selecting all possible pairs
        // real* W = weights.data();
        // for(int p=0; p<N; ++p){
        //     for(int q=0; q<N; ++q){
        //         f_pq = kernel_func(cdf[p]-cdf[q]);
        //         dvar += f_pq * W[p]*W[q]/(W[p]+W[q])*pow(phis[p]-phis[q],2) / N;
        //     }
        // }
        // coeff = variance(phis, weights) / dvar;
        // cout << "coeff" << coeff << endl;

        // set number of mixing pairs
        real num_mix, pmix, alpha, phi_pq;
        num_mix = 3./2. * coeff * N * Omega_phi * dt;
        nmix = floor( num_mix + 1);
        pmix = num_mix / nmix;

        // mixing
        for(int imix=0; imix<nmix; ++imix){
            // select p: marginal prob ~ Wavr + W[p]
            do {
                p = floor(randr()*N);
            } while (W[p] + Wavr <= randr() * (Wmax + Wavr));
            // select q: conditional prob ~ W[p] + W[q]
            do {
                q = floor(randr()*N);
            } while (W[q] + W[p] <= randr() * (Wmax + W[p]));
            if (p==q) continue;
            
            f_pq = kernel_func(cdf[p]-cdf[q]);
            if (randr()<=f_pq && randr()<=pmix){
                alpha = randr();
                phi_pq = (phis[p]*W[p] + phis[q]*W[q])/(W[p]+W[q]);
                phis[p] = phis[p] - alpha * (phis[p] - phi_pq);
                phis[q] = phis[q] - alpha * (phis[q] - phi_pq);
            }
        }
        toc();
    }
};
