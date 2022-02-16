#include "MixingModels.hpp"

// gaussian PDF
real pdfFunc::eval(real x) const {
    return 1./sqrt(2*PI) * exp(-pow(x,2)/2);
}

int main()
{
    /*************************
     mixing settings
     *************************/
    int N = 4000;              // number of samples
    int MAX_STEPS = 1000;       // maximum steps
    real Omega_phi = 2;         // mixing frequency
    real dt = 4e-3;             // timestep
    string PATH = "./data/";    // data directory
    string filename;

    /*************************
     generate samples
     *************************/
    vector<real> particles, weights;

    // // Case F2c of DNS PoF_1988
    // string casename = "PoF_1988_F2c";
    // real max_j = 5;
    // real var_arr[5] = {0.92, 0.80, 0.54, 0.35, 0.28};
    // string PDF_file = "./data/inert/PoF_DNS_1988_F2c_0.92.txt";

    // Case Fig9b of DNS PoF_1996
    string casename = "PoF_1996_Fig9b";
    real max_j = 6;
    real var_arr[6] = {1.0, 0.8, 0.6, 0.5, 0.4, 0.3};
    string PDF_file = "./data/inert/Juneja_1996_PoF_Phi2_Var_1.0.txt";

    /** sampling with given PDF points (X, Y) **/
    vector<real> phi, P;
    readXY(PDF_file, phi, P);
    acceptSampling(phi, P, N, particles, weights, "uniform");

    // /** sampling with PDF function f(X) **/
    // real xlim[2] = {-4, 4};
    // acceptSampling(pdfFunc(xlim), xlim, N, particles, weights, "uniform");

    /*************************
     set mixing models
     *************************/
    vector<MixingModel*> MMs;
    MMs.push_back(new IEM(particles, weights));
    MMs.push_back(new MC(particles, weights));
    MMs.push_back(new EMST(particles, weights));
    MMs.push_back(new KerM(particles, weights, 0.25));

    cout << std::flush;
    /*************************
     mixing for each model
     *************************/
    for(int k=0; k<MMs.size(); ++k){
        MixingModel* mm = MMs[k];
        cout << "Start Mixing " << mm->name << endl;

        int j = 0;
        real std_j, std_0 = sqrt(variance(mm->phis, mm->weights));
        for(int i=0; i<MAX_STEPS; ++i){
            mm->update(Omega_phi, dt);

            std_j = sqrt(variance(mm->phis, mm->weights));
            if (std_j/std_0*var_arr[0] <= var_arr[j]){
                filename = PATH + casename + "_" + mm->name + "_"
                           + to_string(var_arr[j]) + ".txt";
                printf("i = %4d, j = %d\r\n", i, j);
                cout << std::flush;
                saveXY(filename, mm->phis, mm->weights);
                j += 1;
            }
            if(j==max_j) break;
        }
        cout << "Time Cost = " << mm->timecost << endl;
        filename = PATH + casename + "_" + mm->name + "_costs.txt";
        saveCost(filename, N, mm->timecost);
    }
    return 0;
}
