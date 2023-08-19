#include "itensor/all.h"
#include "tdvp.h"
#include "TFIsingLR.h"

using namespace itensor;

//magnetic field vector
std::vector<double> hvector(int, double, double, double, double, double, double);
// calculate local energy density for long range model using innner(psi, hi, psi) 
// with truncation error of the MPO of 1E-6
std::vector<double> calculateLocalEnergy(MPS, std::vector<MPO>);
//calculate Von Neumann entanglement entropy
Real vonNeumannS(MPS, int);
//calculate spin-spin correlator
double spinspin(int,int,MPS,SiteSet);

int main(int argc, char *argv[]){
    std::clock_t tStart = std::clock();

    if(argc < 3){ 
        printfln("Usage: %s input_file",argv[0]); 
        return 0; 
        }
    auto inputSim = InputGroup(argv[1],"input");

    //get simulation parameters
    auto N = inputSim.getInt("N");
    auto v = inputSim.getReal("v", 2.);
    auto h = inputSim.getReal("h", 1.);
    auto hc = inputSim.getReal("hc",1.);
    auto tau = inputSim.getReal("tau", 0.4);
    auto truncE = inputSim.getReal("truncE", 1E-10);
    auto maxB = inputSim.getInt("maxDim", 256);
    auto alpha = inputSim.getReal("alpha",6.);
    auto k = inputSim.getInt("k",3);
    auto dt = inputSim.getReal("dt", 0.1);
    auto tanhshift = inputSim.getReal("tanhshift", 2.);
    
    //get ak, bk
    auto inputExp = InputGroup(argv[2],"input");
    std::vector<double> ak(k), bk(k);
    for (int i = 0; i < k; i++){
        char as[16], bs[16];
        int na = sprintf (as, "a%d", i+1);
        int nb = sprintf (bs, "b%d", i+1);
        ak[i] = inputExp.getReal(as);
        bk[i] = inputExp.getReal(bs);
    }

    // We will write into a file with the time-evolved energy density at all times.
    char schar2[128];
    char schar3[128];
    int n2 = std::sprintf(schar2,"N_%d_v_%0.2f_h_%0.2f_tau_%0.2f_maxDim_%d_alpha_%0.1f_k_%d_tfiLRSuper_En.dat"
                                    ,N,v,h,tau,maxB,alpha,k);
    int n3 = std::sprintf(schar3,"N_%d_v_%0.2f_h_%0.2f_tau_%0.2f_maxDim_%d_alpha_%0.1f_k_%d_tfiLRSuper_SSC.dat"
                                    ,N,v,h,tau,maxB,alpha,k);

    std::string s2(schar2), s3(schar3);
    std::ofstream enerfile, sscfile;
    enerfile.open(s2); // opens the file
    if( !enerfile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    sscfile.open(s3); // opens the file
    if( !sscfile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    enerfile << "time" << " " << "energy" << " " << "SvN" << " " << "bondDim" << " " << "localEnergy" << " " << std::endl;
    sscfile << "time" << " " << "sxsx" << " " << std::endl;
    
    auto sites = SpinHalf(N,{"ConserveQNs=",false,"ConserveParity=",true});
    auto state = InitState(sites);
    for(int i = 1; i <= N; i++){
        state.set(i,"Up");
    }
    auto initState = MPS(state);
    PrintData(totalQN(initState));
    
    // Create the LR Hamiltonian
    std::vector<double> hc_vec(N, hc);
    MPO Hfinal = TFIsingLR(sites, ak, bk, hc_vec);
    
    //sweeps
    auto sweeps = Sweeps(6); //number of sweeps is 5
    sweeps.maxdim() = 10,20,50,100,maxB; //gradually increase states kept
    sweeps.cutoff() = truncE; //desired truncation error
    
    // Create the Local Energy Density Tensors
    std::vector<double> localEnergy(N-1);
    std::vector<MPO> LED(N-1);
    for (int b = 1; b < N; b++){
        auto ampo = AutoMPO(sites);
        ampo += -4.0, "Sx", b, "Sx", b+1;
        ampo += -hc , "Sz", b, "Id", b+1;
        ampo += -hc , "Id", b, "Sz", b+1;
        for( int r = 2; r <= N; r+=2){ // r even
            if( 1 <= b-r/2 && b+r/2+1 <= N ){ // open boundary conditions
                double J = 0.; // get interaction strength
                for (int j = 1; j <= k; j++){
                    J += ak[j-1]*pow( bk[j-1], double(r-1) );
                }
                ampo += -2.0*J, "Sx", b-r/2,   "Sx", b+r/2;
                ampo += -2.0*J, "Sx", b-r/2+1, "Sx", b+r/2+1;
            }
        }
        for( int r = 3; r <= N; r+=2){ // r odd
            if( 1 <= b-(r-1)/2 && b+(r+1)/2 <= N ){ // open boundary conditions
                double J = 0.; // get interaction strength
                for (int j = 1; j <= k; j++){
                    J += ak[j-1]*pow( bk[j-1], double(r-1) );
                }
                ampo += -4.0*J, "Sx", b-(r-1)/2, "Sx", b+(r+1)/2;
            }
        }
        
        // choose cutoff depending on the power of alpha
        // if (N-1)^-alpha > 1E-6, choose N^-alpha as the cutoff
        LED[b-1] = toMPO(ampo, {"Cutoff", 1E-6});
        // printfln("max link dimension of MPO at bond %d = %d", b, maxLinkDim(LED[b-1]));
    }

    // Create the SxSx vector
    std::vector<double> sxsxcorr(N);

    //magnetic field vector
    std::vector<double> hvals = hvector(N, 0.0, h, hc, v, tau, tanhshift);
    //either be able to create an MPO with a vector of h or locally update the MPO. Using function?
    MPO H = TFIsingLR(sites, ak, bk, hvals);
    
    // Find Initial Ground State
    auto [en,psi] = dmrg(H,initState,sweeps,{"Silent=",true});
    en = inner(psi, Hfinal, psi);
    //calculate entanglement
    auto SvN = vonNeumannS(psi, N/2);
    //get bond dimensions
    IndexSet bonds = linkInds(psi); 
    //calculate local energy <psi|Hf(x)|psi>
    localEnergy = calculateLocalEnergy(psi, LED);
    //calculate spin-spin correlation
    for (int b = 1; b <= N; b++){
        sxsxcorr[b-1] = spinspin(N/2,b,psi,sites);
    }

    //store variables to energy file
    enerfile << 0.0 << " " << en << " " << SvN << " ";
    for (int j=0; j<N-1; j++){
        enerfile << dim(bonds[j]) << " ";
    }
    for (int j = 0; j < N-1; j++){
        enerfile << localEnergy[j] << " ";
    }
    enerfile << std::endl;
    //store variables to spin spin correlation file
    sscfile << 0.0 << " ";
    for (int j = 0; j < N; j++){
        sscfile << sxsxcorr[j] << " ";
    }
    sscfile << std::endl;

    // time evolution parameters. Get time accuracy of 1E-4
    double tval = 0.0;
    double delta1 =  0.414490771794376*dt;
    double delta2 = -0.657963087177503*dt;
    double finalTime = 0.5*double(N) + 0.5*double(N)/v + 2.0*tau*tanhshift; // 0.1*N/c + 0.5*N/v + 2*tau*shift
    int nt = int(finalTime/dt);
    int linkCheck = int(log2( double(maxB)) );
    int numCenter = 2; // start with two-site tdvp

    // 4th order TDVP parameters
    auto sweeps1 = Sweeps(2); //two forward time steps of delta1
    sweeps1.maxdim() = maxB;
    sweeps1.cutoff() = truncE;
    sweeps1.niter() = 10;
    auto sweeps2 = Sweeps(1); //one backward time step of delta2
    sweeps2.maxdim() = maxB;
    sweeps2.cutoff() = truncE;
    sweeps2.niter() = 10;

    printfln("t = %0.1f, energy = %0.2f, SvN = %0.3f", tval, en, SvN);

    println("Starting 4th order two-site TDVP");

    ////////////////////
    // TIME EVOLUTION //
    ////////////////////
    for (int n = 1; n <= nt; n++){

        tval += dt;
        
        //get new hvals
        std::vector<double> hvals = hvector(N, tval, h, hc, v, tau, tanhshift);
        //either create new MPO every time or update existing MPO. I think updating is better. Could be done in function
        H = TFIsingLR(sites, ak, bk, hvals);

        // TDVP sweep
        tdvp(psi, H, -Cplx_i*delta1, sweeps1,{"Silent",true,"NumCenter",numCenter});
        tdvp(psi, H, -Cplx_i*delta2, sweeps2,{"Silent",true,"NumCenter",numCenter});
        tdvp(psi, H, -Cplx_i*delta1, sweeps1,{"Silent",true,"NumCenter",numCenter});        
        
        // calculate energy <psi|Hf|psi>
        en = innerC(psi, Hfinal, psi).real();
        //calculate entanglement entropy
        SvN = vonNeumannS(psi, N/2);
        //calculate local energy <psi|Hf(x)|psi>
        localEnergy = calculateLocalEnergy(psi, LED);

        enerfile << tval << " " << en << " " << SvN << " ";
        auto bonds = linkInds(psi); //get bond dimensions
        for (int j = 0; j < N-1; j++){
            enerfile << dim(bonds[j]) << " ";
        }
        for (int j = 0; j < N-1; j++){
            enerfile << localEnergy[j] << " ";
        }
        enerfile << std::endl;

        printfln("t = %0.1f, energy = %0.2f, SvN = %0.3f", tval, en, SvN);

        //calculate spin-spin correlation
        for (int b = 1; b <= N; b++){
            sxsxcorr[b-1] = spinspin(N/2,b,psi,sites);
        }

        //store variables to spin spin correlation file
        sscfile << tval << " ";
        for (int j = 0; j < N; j++){
            sscfile << sxsxcorr[j] << " ";
        }
        sscfile << std::endl;


        // check if bondDim is maxed out
        if( numCenter > 1 && dim(bonds[linkCheck-1]) >= maxB){
            printfln("link %d has bond dimension %d", linkCheck, dim(bonds[linkCheck-1]));
            printfln("Switching to 4th order 1-site TDVP");
            numCenter = 1;
        }

    }// for n
    
    enerfile.close();
    sscfile.close();

    print(" END PROGRAM. TIME TAKEN :");
    printf("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
}

std::vector<double> hvector(int N, double tval, double h, double hc, double v, double tau, double tanhshift)
    {
    std::vector<double> hvals(N);
    for (int b = 1; b <= N; b++){
        double f = abs(double(b-N/2)-0.5)/v - tval;
        hvals[b-1] = hc + h*(0.5 + 0.5*tanh( f/tau + tanhshift ));
    }
    return hvals;
}

std::vector<double> calculateLocalEnergy(MPS psi, std::vector<MPO> LED){

    int N = psi.length();
    std::vector<double> energy(N-1);

    for( int b = 1; b < N; b++){
        
        energy[b-1] = innerC(psi, LED[b-1], psi).real();

    }

    return energy;
}

//calculate entanglement
Real vonNeumannS(MPS psi, int b){
    Real SvN = 0.;

    //choose orthogonality center and perform svd
    psi.position(b);
    auto l = leftLinkIndex(psi,b);
    auto s = siteIndex(psi,b);
    auto [U,S,V] = svd(psi(b),{l,s});
    auto u = commonIndex(U,S);

    //Apply von Neumann formula
    //to the squares of the singular values
    for(auto n : range1(dim(u))){
        auto Sn = elt(S,n,n);
        auto p = sqr(Sn);
        if(p > 1E-12) SvN += -p*log(p);
    }
    return SvN;

}//vonNeumannS

//calculate spin-spin correlator
double spinspin(int center, int b, MPS psi, SiteSet sites){
    
    double corrX;

    psi.position(b);
    if(b>center){ //bring site b next to the center from right
        for(int n=b-1; n>center; n--){
            auto g = BondGate(sites,n,n+1);
            auto AA = psi(n)*psi(n+1)*g.gate(); //contract over bond n
            AA.replaceTags("Site,1","Site,0");
            psi.svdBond(g.i1(), AA, Fromright); //svd from the right
            psi.position(g.i1()); //move orthogonality center to the left 
        }
        auto ket = psi(center)*psi(center+1);
        auto SxSx = 4.0*sites.op("Sx",center)*sites.op("Sx",center+1);
        corrX = eltC( dag(prime(ket,"Site")) * SxSx * ket).real();
    }
    else if(b<center){ //bring site b next to the center from left
        for(int n=b; n<center-1; n++){
          auto g = BondGate(sites,n,n+1);
          auto AA = psi(n)*psi(n+1)*g.gate(); //contract over bond n
          AA.replaceTags("Site,1","Site,0");
          psi.svdBond(g.i1(), AA, Fromleft); //svd from the right
          psi.position(g.i2()); //move orthogonality center to the right 
        }
        auto ket = psi(center-1)*psi(center);
        auto SxSx = 4.0*sites.op("Sx",center-1)*sites.op("Sx",center);
        corrX = eltC( dag(prime(ket,"Site")) * SxSx * ket).real();
    }
    else{
        corrX = 1.;
    }

    return corrX;

}//SxSx
