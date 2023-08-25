//
// Copyright [2023] [Simon Bernier]
//
#include "itensor/all.h"
#include "TFIsingLR.h"

using namespace itensor;

//magnetic field vector
std::vector<double> hvector(int, double, double, double);
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
    auto hc = inputSim.getReal("hc",1.);
    auto truncE = inputSim.getReal("truncE", 1E-10);
    auto maxB = inputSim.getInt("maxDim", 256);
    auto alpha = inputSim.getReal("alpha",6.);
    auto k = inputSim.getInt("k",3);
    auto h = inputSim.getReal("h", 4.);
    auto R = inputSim.getReal("R", 1.);
    
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
    char schar[128];
    int n1 = std::sprintf(schar,"N_%d_h_%0.2f_R_%0.1f_alpha_%0.1f_k_%d_tfiLRstatic.dat",N,h,R,alpha,k);

    std::string s1(schar);
    std::ofstream datafile;
    datafile.open(s1); // opens the file
    if( !datafile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    datafile << "energy" << " " << "var" << " " << "SvN" << " " << "localEnergy" << " " << "sxsx" << " " << std::endl;
    
    auto sites = SpinHalf(N,{"ConserveQNs=",false,"ConserveParity=",true});
    auto state = InitState(sites);
    for(int i = 1; i <= N; i++){
        state.set(i,"Up");
    }
    auto initState = MPS(state);
    PrintData(totalQN(initState));
    
    // Create the LR Hamiltonian
    std::vector<double> hc_vec = hvector(N, h, hc, R);
    MPO H = TFIsingLR(sites, ak, bk, hc_vec);
    
    //sweeps
    auto sweeps1 = Sweeps(15); //number of sweeps is 5
    sweeps1.maxdim() = 10,20,50,100,maxB; //gradually increase states kept
    sweeps1.cutoff() = truncE; //desired truncation error
    sweeps1.noise()  = 1E-7, 1E-8, 0, 0, 1E-7, 1E-8, 0;

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
        LED[b-1] = toMPO(ampo,  {"Cutoff", 1E-6 });
    }

    // Create the SxSx vector
    std::vector<double> sxsxcorr(N);
   
    // Find Initial Ground State
    auto [energy,psi] = dmrg(H,initState,sweeps1,{"Silent=",true});
    auto var = inner(H,psi,H,psi) - energy*energy;
    printfln("var = %0.3g", var);

    //calculate entanglement
    auto SvN = vonNeumannS(psi, N/2+1);
    //calculate local energy <psi|Hf(x)|psi>
    localEnergy = calculateLocalEnergy(psi, LED);
    //calculate spin-spin correlation
    for (int b = 1; b <= N; b++){
        sxsxcorr[b-1] = spinspin(N/2,b,psi,sites);
    }

    //store variables to energy file
    datafile << energy << " " << var << " " << SvN << " ";
    for (int j = 0; j < N-1; j++){
        datafile << localEnergy[j] << " ";
    }
    //store variables to spin spin correlation file
    for (int j = 0; j < N; j++){
        datafile << sxsxcorr[j] << " ";
    }
    datafile << std::endl;

    datafile.close();

    printfln("energy = %0.3f, SvN = %0.3f, maxDim = %d", energy, SvN, maxLinkDim(psi));

    print(" END PROGRAM. TIME TAKEN :");
    printf("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
}

// make h-vector
std::vector<double> hvector(int N, double h, double hc, double R){

    std::vector<double> hvals(N);

    if (R==0.0){
        for (int b = 1; b <= N/4; b++){
            hvals[b-1] = hc+h;
        }
        for (int b = N/4+1; b <= 3*N/4; b++){
            hvals[b-1] = hc;
        }
        for (int b = 3*N/4+1; b <= N; b++){
            hvals[b-1] = hc+h;
        }
    }
    else{
        for (int b = 1; b <= N; b++){
            double f = abs(double(b-N/2)-0.5)-N/4;
            hvals[b-1] = hc + h*(0.5 + 0.5*tanh( f/R ));
        }
    }
    
    return hvals;
}

// calculate energy using MPO
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
        corrX = elt( dag(prime(ket,"Site")) * SxSx * ket);
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
        corrX = elt( dag(prime(ket,"Site")) * SxSx * ket);
    }
    else{
        corrX = 1.;
    }

    return corrX;

}//SxSx
