//
// Copyright [2023] [Simon Bernier]
//
#include "itensor/all.h"
#include "TFIsingLR.h"

using namespace itensor;

int main(int argc, char *argv[]){
    std::clock_t tStart = std::clock();

    if(argc < 2){ 
        printfln("Usage: %s input_file",argv[0]); 
        return 0; 
    }
    auto input = InputGroup(argv[1],"input");

    auto k = input.getInt("k");
    auto alpha = input.getReal("alpha");
    std::vector<double> ak(k), bk(k);
    for (int i = 0; i < k; i++){
        char as[16], bs[16];
        int na = sprintf (as, "a%d", i+1);
        int nb = sprintf (bs, "b%d", i+1);
        ak[i] = input.getReal(as);
        bk[i] = input.getReal(bs);
    }

    int N=16;
    double h=0.;
    if(argc > 2)
        N = std::stoi(argv[2]);
    if(argc > 3)
        h = std::stod(argv[3]);
    
    printfln("N = %d, h = %0.2f", N, h);

    char schar1[128];
    int n1 = std::sprintf(schar1,"N_%d_h_%0.2f_alpha_%0.1f_k_%d_tfi_gap.dat", N, h, alpha, k);

    std::string s1(schar1);
    std::ofstream enerfile;
    enerfile.open(s1); // opens the file
    if( !enerfile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    enerfile << "e0" << " " << "var0" << " " << "e1" << " "  << "var1" << " " << std::endl;
    
    auto sites = SpinHalf(N,{"ConserveQNs=",false,"ConserveParity=",true});

    auto state = InitState(sites);
    for(int i = 1; i <= N; i++){
        state.set(i,"Up");
    }
    auto initState = MPS(state);
    PrintData(totalQN(initState));
    
    // Create the Heisenberg Hamiltonian
    MPO H = TFIsingLR(sites, ak, bk, {"h",h});
    //PrintData(H);
    
    //sweeps
    auto sweeps = Sweeps(15); //number of sweeps is 6
    sweeps.maxdim() = 10,20,50,100,100,200; //gradually increase states kept
    sweeps.cutoff() = 1E-10; //desired truncation error
    sweeps.noise() = 1E-7,1E-8,0,1E-8,0,1E-8,0;
    
    // Find Initial Ground State
    auto [en0,psi0] = dmrg(H,initState,sweeps,{"Silent=",true});
    auto var0 = inner(H, psi0, H, psi0) - en0*en0;
    printfln("E0 = %0.8f, var0 = %0.8g", en0, var0);

    // find first excited state
    auto wfs = std::vector<MPS>(1);
    wfs.at(0) = psi0;
    auto [en1,psi1] = dmrg(H,wfs,initState,sweeps,{"Silent=",true,"Weight=",20.0});
    auto var1 = inner(H,psi1,H,psi1)-en1*en1;
    printfln("E1 = %0.8f, var1 = %0.8g", en1, var1);

    printfln("\n gap = %0.3g", en1-en0);
    printfln("overlap = %0.3g", inner(psi1,psi0));

    //store variables to energy file
    enerfile << en0 << " " << var0 << " " << en1 << " " << var1 << " " << std::endl;
    enerfile.close();

    print(" END OF PROGRAM. ");
    printf("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
}
