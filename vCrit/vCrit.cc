#include "itensor/all.h"
#include "tdvp.h"
#include "TFIsingLR.h"

using namespace itensor;

//calculate Von Neumann entanglement entropy
Real vonNeumannS(MPS, int);

int main(int argc, char *argv[]){
    std::clock_t tStart = std::clock();

    if(argc < 2){ 
        printfln("Usage: %s input_file",argv[0]); 
        return 0; 
    }
    auto input = InputGroup(argv[1],"input");

    auto N = input.getInt("N", 16);
    auto h = input.getReal("h", 4.0);
    auto k = input.getInt("k");
    auto alpha = input.getReal("alpha");
    auto truncE = input.getReal("truncE", 1E-10);
    auto maxB = input.getInt("maxDim", 128);
    std::vector<double> ak(k), bk(k);
    for (int i = 0; i < k; i++){
        char as[16], bs[16];
        int na = sprintf (as, "a%d", i+1);
        int nb = sprintf (bs, "b%d", i+1);
        ak[i] = input.getReal(as);
        bk[i] = input.getReal(bs);
    }
    
    // We will write into a file with the time-evolved energy density at all times.
    char schar1[128];
    int n1 = std::sprintf(schar1,"N_%d_h_%0.2f_alpha_%0.1f_k_%d_TFIvCrit.dat", N,h,alpha,k);

    std::string s1(schar1);
    std::ofstream dataFile;
    dataFile.open(s1); // opens the file
    if( !dataFile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    
    //make header for t=0 calculations
    dataFile << "t=0" << " " << "enPsi" << " " << "enPhi" << " " << "svn(x)" << " " << "bondDim" << " " << std::endl;

    // make vectors to store SvN
    std::vector<Real> svn(N-1,0.);
    
    auto sites = SpinHalf(N,{"ConserveQNs=",false,"ConserveParity=",true});
    auto state = InitState(sites);
    for(int i = 1; i <= N; i++){
        state.set(i,"Up");
    }
    auto initState = MPS(state);
    PrintData(totalQN(initState));
    
    // Create the LR Hamiltonian
    MPO H = TFIsingLR(sites, ak, bk, {"h",h});
    
    //sweeps
    auto sweeps = Sweeps(10); //number of sweeps
    sweeps.maxdim() = 10,20,50,100,maxB; //gradually increase states kept
    sweeps.cutoff() = truncE; //desired truncation error
    sweeps.noise() = 1E-7, 1E-8, 0;
    
    // Find Initial Ground State
    auto [en_psi,psi] = dmrg(H,initState,sweeps,{"Silent=",true});

    // make |phi> = sigma_z|psi>
    int loc = N/2; //centered in the middle of the chain
    psi.position(loc);
    auto newA = 2.0*sites.op("Sz",loc)*psi(loc);
    newA.noPrime();
    psi.set(loc, newA);
    psi.orthogonalize({"Cutoff=",truncE,"MaxDim=",maxB});
    auto en_phi = innerC(psi,H,psi).real(); //energy after disturbing the ground state

    for(auto b : range1(N-1)){
        svn[b-1] = vonNeumannS(psi, b);
    }
    auto bonds = linkInds(psi); //get bond dimensions

    printfln("\nt=0; phi energy = %0.3f, max link dim is %d", en_phi, maxLinkDim(psi));
    // store to file
    dataFile << 0.0 << " " << en_psi << " " << maxLinkDim(psi) << " " << en_phi << " " ;
    for(int j = 0; j<N-1; j++){ //save svn
        dataFile << svn[j] << " ";
    }
    for(int j = 0; j<N-1; j++){ //save bond dim
        dataFile << dim(bonds[j]) << " ";
    }
    dataFile << std::endl;

    //make header for t>0 calculations
    dataFile << "t" << " " << "enPhi" << " " << "svn(x)" << " " << "bondDim" << " " << std::endl;
 
    // time evolution parameters.
    double tval = 0., dt = 0.1;
    double delta1 =  0.414490771794376*dt; // for 4th order TDVP
    double delta2 = -0.657963087177503*dt;
    double finalTime = 32.;
    int nt=int(finalTime/dt);
    int linkCheck = int(log2( double(maxB)) );
    int numCenter = 2; // start with two-site tdvp

    // 4th order TDVP parameters
    auto sweeps1 = Sweeps(2); //two forward time steps of delta1
    sweeps1.maxdim() = maxB;
    sweeps1.cutoff() = truncE;
    sweeps1.niter() = 15;
    auto sweeps2 = Sweeps(1); //one backward time step of delta2
    sweeps2.maxdim() = maxB;
    sweeps2.cutoff() = truncE;
    sweeps2.niter() = 15;

    println("Starting 4th order two-site TDVP");
    
    ////////////////////
    // TIME EVOLUTION //
    ////////////////////
    for (int n = 1; n <= nt; ++n){
        tval += dt;

        // TDVP sweep
        tdvp(psi, H, -Cplx_i*delta1, sweeps1, {"Silent",true,"NumCenter",numCenter});
        tdvp(psi, H, -Cplx_i*delta2, sweeps2, {"Silent",true,"NumCenter",numCenter});
        en_phi = tdvp(psi, H,- Cplx_i*delta1, sweeps1, {"Silent",true,"NumCenter",numCenter});        
        
        // calculate svn
        for(auto b : range1(N-1)){
            svn[b-1] = vonNeumannS(psi, b);
        }
        bonds = linkInds(psi); //get bond dimensions

        //write to file
        dataFile << tval << " " << en_psi << " " << maxLinkDim(psi) << " ";
        for(int j = 0; j<N-1; j++){ //save svn
            dataFile << svn[j] << " ";
        }
        for(int j = 0; j<N-1; j++){ //save bond dim
            dataFile << dim(bonds[j]) << " ";
        }
        dataFile << std::endl;

        printfln("\nIteration %d, time = %0.2f; phi energy = %0.3f, max link dim is %d",n,tval,en_phi,maxLinkDim(psi));

        // check if bondDim is maxed out
        if( numCenter > 1 && dim(bonds[linkCheck-1]) >= maxB){
            printfln("\nlink %d has bond dimension %d", linkCheck, dim(bonds[linkCheck-1]));
            printfln("Switching to 4th order 1-site TDVP");
            numCenter = 1;
        }

    }// for n
    
    dataFile.close();

    print(" END OF PROGRAM. ");
    printf("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
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