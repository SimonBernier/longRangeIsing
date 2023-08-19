//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HAMS_TFISINGLR_H
#define __ITENSOR_HAMS_TFISINGLR_H

#include "itensor/mps/mpo.h"

namespace itensor {

class TFIsingLR
    {
    public:

    TFIsingLR(SiteSet const& sites, 
              std::vector<double> const& ak,
              std::vector<double> const& bk,
              std::vector<double> const& h,
              Args const& args = Args::global());

    operator MPO() { init_(); return H; }

    private:

    //////////////////
    //
    // Data Members

    SiteSet const& sites_;
    std::vector<double> const& ak_, //prefactor
                               bk_, //exponential decay
                               h_; // transverse field
    int N_, // system size
        k_; // number of exponentials
    Real J_; // hopping
    bool initted_, // no idea
         infinite_; // no idea
    MPO H;

    //
    //////////////////

    void 
    init_();

    }; //class TFIsingLR

inline TFIsingLR::
TFIsingLR(SiteSet const& sites, 
          std::vector<double> const& ak,
          std::vector<double> const& bk,
          std::vector<double> const& h,
          Args const& args)
  : sites_(sites),
    ak_(ak),
    bk_(bk),
    h_(h),
    initted_(false)
    { 
    N_ = sites_.length();
    k_ = ak_.size();
    J_ = args.getReal("J",1.);
    infinite_ = args.getBool("Infinite",false);
    }

void inline TFIsingLR::
init_()
    {
    if(initted_) return;

    H = MPO(sites_);

    std::vector<Index> links(N_+1);

    for(int l = 0; l <= N_; ++l) 
        {
        auto ts = format("Link,l=%d",l);
        links.at(l) = Index(QN({"Parity",1,2}),2,
                            QN({"Parity",0,2}),k_,
                            Out,
                            ts);
        }

    Index const& last = (infinite_ ? links.at(0) : links.at(N_));

    for(int n = 1; n <= N_; ++n)
        {
        auto& W = H.ref(n);
        auto row = dag(links.at(n-1));
        auto col = (n==N_ ? last : links.at(n));

        W = ITensor(dag(sites_(n)),prime(sites_(n)),row,col);

        W += sites_.op("Id",n) * setElt(row(1)) * setElt(col(1)); //ending state
        W += sites_.op("Id",n) * setElt(row(2)) * setElt(col(2)); //starting state

        for( int j = 0; j < k_; j++){ //exponential decay
            W += sites_.op("Id",n) * setElt(row(3+j)) * setElt(col(3+j)) * bk_[j]; // bk Id
        }

        // problem within this for-loop
        for( int j = 0; j<k_; j++){ //long-range interactions 
            W += sites_.op("Sx",n) * setElt(row(3+j)) * setElt(col(1)) * 2.;                // Sx
            W += sites_.op("Sx",n) * setElt(row(2)) * setElt(col(3+j)) * -2. * J_ * ak_[j]; // ak Sx
        }

        W += sites_.op("Sz",n) * setElt(row(2)) * setElt(col(1)) * -2.*h_[n-1]; //transverse field

        }
    
    auto LH = setElt(links.at(0)(2));
    auto RH = setElt(dag(last)(1));

    if(not infinite_)
        {
        //Multiply first and last
        //MPO tensor by edge vectors
        H.ref(1) *= LH;
        H.ref(N_) *= RH;
        }
    else
        {
        //Store edge vectors just before
        //and after first and last sites
        H.ref(0) = LH;
        H.ref(N_+1) = RH;
        }

    initted_ = true;
    }

} //namespace itensor

#endif
