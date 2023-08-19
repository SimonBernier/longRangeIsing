#include "itensor/all.h"
#include "TFIsingLR.h"

using namespace itensor;

// calculate local energy density for long range model using innner(psi, hi, psi) 
// with truncation error of the MPO of 1E-6
std::vector<double> calculateLocalEnergy(MPS, std::vector<MPO>);
/*
//calculates local energy density
std::vector<double> calculateLocalEnergy(int, MPS, double, std::vector<ITensor>, std::vector<double>, SiteSet);
// calculate long range energy density
double calculateLRenergy(int, int, int, MPS, std::vector<ITensor>, std::vector<double>, SiteSet);
*/
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
    int n2 = std::sprintf(schar2,"N_%d_h_%0.2f_alpha_%0.1f_k_%d_cutoff1E-6_tfiLRcrit_En.dat",N,hc,alpha,k);
    int n3 = std::sprintf(schar3,"N_%d_h_%0.2f_alpha_%0.1f_k_%d_cutoff1E-6_tfiLRcrit_SC.dat",N,hc,alpha,k);

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
    enerfile << "energy" << " " << "var" << " " << "SvN" << " " << "bondDim" << " " << "localEnergy" << " " << std::endl;
    sscfile << "sxsx" << " " << std::endl;
    
    auto sites = SpinHalf(N,{"ConserveQNs=",false,"ConserveParity=",true});
    auto state = InitState(sites);
    for(int i = 1; i <= N; i++){
        state.set(i,"Up");
    }
    auto initState = MPS(state);
    PrintData(totalQN(initState));
    
    // Create the LR Hamiltonian
    std::vector<double> hc_vec(N, hc);
    MPO H = TFIsingLR(sites, ak, bk, hc_vec);
    
    //sweeps
    auto sweeps = Sweeps(10); //number of sweeps is 5
    sweeps.maxdim() = 10,20,50,100,maxB; //gradually increase states kept
    sweeps.cutoff() = truncE; //desired truncation error
    sweeps.noise()  = 1E-7, 1E-8, 0;

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
        //printfln("max link dimension of MPO at bond %d = %d", b, maxLinkDim(LED[b-1]));
    }

    // Create the SxSx vector
    std::vector<double> sxsxcorr(N);
   
    // Find Initial Ground State
    auto [energy,psi] = dmrg(H,initState,sweeps,{"Silent=",true});
    auto var = inner(H,psi,H,psi) - energy*energy;
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
    enerfile << energy << " " << var << " " << SvN << " ";
    for (int j=0; j<N-1; j++){
        enerfile << dim(bonds[j]) << " ";
    }
    for (int j = 0; j < N-1; j++){
        enerfile << localEnergy[j] << " ";
    }
    enerfile << std::endl;
    //store variables to spin spin correlation file
    for (int j = 0; j < N; j++){
        sscfile << sxsxcorr[j] << " ";
    }
    sscfile << std::endl;

    enerfile.close();
    sscfile.close();

    printfln("energy = %0.3f, SvN = %0.3f, maxDim = %d", energy, SvN, maxLinkDim(psi));

    print(" END PROGRAM. TIME TAKEN :");
    printf("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
}

std::vector<double> calculateLocalEnergy(MPS psi, std::vector<MPO> LED){

    int N = psi.length();
    std::vector<double> energy(N-1);

    for( int b = 1; b < N; b++){
        
        energy[b-1] = innerC(psi, LED[b-1], psi).real();

    }

    return energy;
}

/*
//calculate local energy density and return a vector of doubles
std::vector<double> calculateLocalEnergy(int N, MPS psi, double hc, std::vector<ITensor> LED, std::vector<double> J, SiteSet sites){
    
    std::vector<double> energy(N-1);

    for (int b=1; b<N; b++){
        psi.position(b);
        auto ket = psi(b)*psi(b+1);
        auto hi = LED[b-1];
        hi += -hc*sites.op("Sz",b)*sites.op("Id",b+1);
        hi += -hc*sites.op("Id",b)*sites.op("Sz",b+1);
        energy[b-1] = eltC( dag(prime(ket,"Site")) * hi * ket).real();
    }

    for(int i = 1; i<int(size(J)); i++){

        std::vector<double> energyLR(N-1-i);

        for(int b=1; b<N-i; b++){
            energyLR[b-1] = calculateLRenergy(N, b, i, psi, LED, J, sites);
        }

        if( i%2 == 0 ){
            for(int b=1; b<N-i; b++){
                energy[b-1+i/2] += energyLR[b-1];
            } //for b
        }
        else{
            for(int b=1; b<N-i-1; b++){
                energy[b-1 + i/2+1] += 0.5 * (energyLR[b-1] + energyLR[b]); // interpolate
            } //for b
        }// if

    }// for i

    return energy;
}//calculateLocalEnergy

double calculateLRenergy(int N, int b, int i, MPS psi, std::vector<ITensor> LED, std::vector<double> J, SiteSet sites){

    double energy=0.;

    if (b+i+1 <= N){

        for (int j=i; j>=1; j--){// switch sites for long-range interaction

            int ind = b+j;
            psi.position(ind);
            auto g = BondGate(sites,ind,ind+1);
            auto AA = psi(ind)*psi(ind+1)*g.gate(); //contract over bond ind
            AA.replaceTags("Site,1","Site,0");
            psi.svdBond(g.i1(), AA, Fromleft); //svd from the left
            psi.position(g.i1()); //restore orthogonality center to the left

        } // for j

        psi.position(b);
        auto ket = psi(b)*psi(b+1);
        energy = J[i]*eltC( dag(prime(ket,"Site")) * LED[b-1] * ket).real();
                    
        for (int j=1; j<=i; j++){// SMART switch sites for next-nearest neighbour interaction

            int ind = b+j;

            psi.position(ind);
            auto g = BondGate(sites,ind,ind+1);
            auto AA = psi(ind)*psi(ind+1)*g.gate(); //contract over bond ind
            AA.replaceTags("Site,1","Site,0");
            psi.svdBond(g.i1(), AA, Fromleft); //svd from the left
            psi.position(g.i2()); //move orthogonality center to the right

        } // for j
    }// if

    return energy;
}//calculateLRenergy
*/

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
