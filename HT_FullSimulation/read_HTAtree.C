#include <iostream>

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

static constexpr Double_t MagneticField = 5.;
static constexpr Double_t K = 0.000299792458*MagneticField;

TH1F h_muon_pt("muon_pt","muon p_{T};p_{T} [GeV]", 100., 0., 100);
TH1F h_muon_phi("muon_phi","muon #phi;#phi [rad]", 100., -3.1415, 3.1415);
TH1F h_muon_theta("muon_theta","muon #theta;#theta [rad]",100., 0., 3.1415);
TH1F h_muon_vz("muon_vz","muon vz;vz [mm]", 100, -10., 10.);

TH2F h2_yx("h2_yx","hit y vs x; x [mm]; y [mm]",500,-1600.,1600.,500., -1600.,1600.);
TH2F h2_rz("h2_rz","hit #rho vs z; z [mm]; #rho [mm]",500,-2500.,2500.,500.,0.,1600.);

TH2F h2_yx_loc("h2_yx_loc","hit y_{loc} vs x_{loc}; x_{loc} [mm]; y_{loc} [mm]",100,-20.,20., 100., -20.,20.);


void read_HTAtree(){

  // --- open the root file
  TFile *input_file = TFile::Open("ntu_muongun_pt3_phi0-30_theta10-170_10k.root");

  // --- retrieve the tree
  TTree *tree = (TTree*) input_file->Get("HTAtree");

  // --- define branch variables and set the branch addresses

  // MC particle
  int   part_pdg;
  float part_px;
  float part_py;
  float part_pz;
  float part_e;
  float part_vx;
  float part_vy;
  float part_vz;
  float part_q;

  tree->SetBranchAddress("part_pdg", &part_pdg);
  tree->SetBranchAddress("part_px",  &part_px);
  tree->SetBranchAddress("part_py",  &part_py);
  tree->SetBranchAddress("part_pz",  &part_pz);
  tree->SetBranchAddress("part_e",   &part_e);
  tree->SetBranchAddress("part_vx",  &part_vx);
  tree->SetBranchAddress("part_vy",  &part_vy);
  tree->SetBranchAddress("part_vz",  &part_vz);
  tree->SetBranchAddress("part_q",   &part_q);

  // tracker hits
  int n_hit;
  std::vector<int> *hit_index = nullptr;
  std::vector<int> *hit_id0 = nullptr;
  std::vector<float> *hit_x = nullptr;
  std::vector<float> *hit_y = nullptr;
  std::vector<float> *hit_z = nullptr;
  std::vector<float> *hit_t = nullptr;
  std::vector<float> *hit_xloc = nullptr;
  std::vector<float> *hit_yloc = nullptr;
 
  tree->SetBranchAddress("n_hit", &n_hit);                                                                                          
  tree->SetBranchAddress("hit_index", &hit_index);
  tree->SetBranchAddress("hit_id0", &hit_id0);
  tree->SetBranchAddress("hit_x", &hit_x);
  tree->SetBranchAddress("hit_y", &hit_y);
  tree->SetBranchAddress("hit_z", &hit_z);
  tree->SetBranchAddress("hit_t", &hit_t);
  tree->SetBranchAddress("hit_xloc", &hit_xloc);
  tree->SetBranchAddress("hit_yloc", &hit_yloc);

  // --- loop over all entries
  Long64_t nEntries = tree->GetEntries();
  for (Long64_t ientry = 0; ientry < nEntries; ++ientry) {

    tree->GetEntry(ientry);
    
    // loop over the MC partiles
    const float part_pt = sqrt(part_px*part_px + part_py*part_py);
    const float part_phi = atan2(part_py, part_px);
    const float part_theta = atan2(part_pt, part_pz);

    h_muon_pt.Fill(part_pt);
    h_muon_phi.Fill(part_phi);
    h_muon_theta.Fill(part_theta);
    h_muon_vz.Fill(part_vz);
    
    cout << part_pt << " " << part_phi << " " << part_theta << " "  << part_vz << endl ;
    cout << n_hit << endl;
      

    // loop over the tracker hits
    for (int ihit=0; ihit<n_hit; ++ihit){

      // CellID encoding: "system:5,side:-2,layer:6,module:11,sensor:8"
      //  system:
      //    VXD barrel = 1
      //    VXD endcap = 2
      //    IT barrel  = 3
      //    IT endcap  = 4
      //    OT barrel  = 5
      //    OT endcap  = 6
      
      const unsigned int system = (unsigned) ( (*hit_id0)[ihit] & 0x1f );
      const int side = (int) ( ((*hit_id0)[ihit] >> 5) & 0x3 );
      const unsigned int layer = (unsigned) ( ((*hit_id0)[ihit] >> 7) & 0x3f );
      
      const float hit_rho = sqrt((*hit_x)[ihit]*(*hit_x)[ihit]+(*hit_y)[ihit]*(*hit_y)[ihit]);

      h2_yx.Fill((*hit_x)[ihit], (*hit_y)[ihit]);
      h2_rz.Fill((*hit_z)[ihit], hit_rho);

      h2_yx_loc.Fill((*hit_xloc)[ihit], (*hit_yloc)[ihit]);

      cout << "    " << ihit << " " <<  system << " " << side << " " << layer << " "
	   << (*hit_x)[ihit] << " " <<  (*hit_y)[ihit] << " " <<  (*hit_z)[ihit] << " " <<  (*hit_t)[ihit] << endl;
      
    } // ihit


    // loop over the tracks
    //for (int itrk=0; itrk<n_trk; ++itrk){
    //
    //  const float trk_pt = K/fabs(trk_curv[itrk]);
    //  const float trk_theta = (TMath::PiOver2()-TMath::ATan(trk_tanl[itrk]));
    //
    //  cout << "TRACK: " <<  trk_pt << " " << trk_phi[itrk] << " " << trk_theta << endl;
    //
    //} // itrk

    
  } // ientry loop
  
  
  // --- close the root file
  input_file->Close();

}
