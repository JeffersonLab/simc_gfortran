#ifndef RECON_HCANA_H
#define RECON_HCANA_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRotation.h"

class recon_hcana {
public:
  // Pass "sidis", "heep", "rho", "delta", "exclusive" for reaction_str
  recon_hcana(TString filename,
		TString reaction_str,
		TString hadron_type="mpi",
		Bool_t Earm_HMS = kTRUE);
  ~recon_hcana();

  // Convert reaction_str to lower cases
  inline void ReadReaction(TString& str){
    for (int i=0;i<str.Length();++i){
        str[i]=tolower(str[i]);
        reaction=str;}
  }

  //Auxiliary Function (from hcana) to calculate Pmx, Pmy, Pmz in the Lab/q-frame correctly
  void GeoToSph(Double_t th_geo, Double_t ph_geo, Double_t& th_sph, Double_t& ph_sph);
  void SetCentralAngles(Double_t th_cent, Double_t ph_cent);
  void TransportToLab(Double_t p, Double_t xptar, Double_t yptar, TVector3& pvect);

  std::string              getString(char x);
  std::vector<std::string> FindString(TString keyword, TString fname);
  std::vector<std::string> split(std::string str, char del=':');
  std::vector<float>       num_split(std::string str);

  // Build paths; .hist under ../../outfiles/, .root under ../../worksim/
  inline void buildFileName(TString filename) {
    InSIMCFilename  = filename; // bare stem
    InSIMCHistname  = "../../outfiles/" + InSIMCFilename + ".hist";
    InSIMCRootname  = "../../worksim/"  + InSIMCFilename + ".root";
    OutSIMCRootname = "../../worksim/recon_hcana_" + InSIMCFilename + ".root";
  }

  // Functions that parse insput, writes output
  void ReadTree(const TString& reaction); // Given reaction type, decides which branches to read/write
  void CommonTree(); // read: common in all reaction type
  void CommonNewTree(); // write: common in all reaction type
  void SRDECommonTree(); // read: common in SidisRhoDeltaExclusive types
  void SRDECommonNewTree(); // write: common in SidisRhoDeltaExclusive types
  void EventLoop(); // Calculations written here. Fill().
  void WriteTree(); // Write to the output file.

  // Input/Output
  TFile *fin  = nullptr;
  TFile *fout = nullptr;
  TTree *tree = nullptr;
  TTree *newTree = nullptr;

  // variable declaration
  Int_t   nentries = 0;
  TString reaction;
  TString InSIMCFilename, InSIMCHistname, InSIMCRootname, OutSIMCRootname;
  Float_t progress = 0.0f;
  Float_t Weight = 0.0f;
  Float_t ErrorVal = -999.0f; // Placeholde value for missing branches
  Bool_t  fEarm_HMS;

  // from .hist (constants per file)
  Int_t    simc_nevents    = 0;      // Ngen
  Double_t simc_normfactor = 0.0;    // normfac
  Double_t Ein             = 0.0;    // MeV. Beam Energy
  Double_t kf0             = 0.0;    // MeV. e-arm momentum
  Double_t Pf0             = 0.0;    // MeV. P-arm momentum
  Double_t e_th=0.0, e_ph=0.0;       // degrees. e-arm central angle
  Double_t h_th=0.0, h_ph=0.0;       // degrees. P-arm central angle

  // detector branches
  Float_t hsdelta, hsyptar, hsxptar, hsytar;
  Float_t hsxfp, hsxpfp, hsyfp, hsypfp, hsdeltai, hsyptari, hsxptari, hsytari;
  Float_t ssdelta, ssyptar, ssxptar, ssytar;
  Float_t ssxfp, ssxpfp, ssyfp, ssypfp, ssdeltai, ssyptari, ssxptari, ssytari;

  // Choose which spectrometer is the electron arm vs hadron arm
  Float_t edelta, exptar, eyptar;  // electron arm optics
  Float_t hdelta, hxptar, hyptar;  // hadron  arm optics

  // keeping physics variables computed in simc-way as var_simc
  // but in the tree it will be just var
  Float_t q_simc=0, nu_simc=0, Q2_simc=0, W_simc=0, epsilon_simc=0;
  Float_t z_simc=0, pt2_simc=0, xbj_simc=0;
  Float_t Em_simc=0, Pm_simc=0, thetapq_simc=0, phipq_simc=0, missmass_simc=0;
  Float_t Pmx_simc=0, Pmy_simc=0, Pmz_simc=0;

  // keeping physics variables computed in hcana-way as var_recon
  Float_t q_recon=0, nu_recon=0, Q2_recon=0, W_recon=0, epsilon_recon=0;
  Float_t z_recon=0, pt2_recon=0, pt_recon=0, ptx_recon=0, pty_recon=0;
  Float_t Em_recon=0, Pm_recon=0, thetapq_recon=0, phipq_recon=0, missmass_recon=0;
  Float_t Pmx_recon=0, Pmy_recon=0, Pmz_recon=0, xbj_recon=0;
  Float_t fWeight=0; // Weight * normfac / Ngen

  // Other variables
  Float_t epscm=0, fry=0;
  Float_t mmnuc=0, phad=0, t=0, u=0, s=0;
  Float_t pmpar=0, pmper=0, pmoop=0, pfermi=0, siglab=0;
  Float_t sigcm=0, decdist=0, Mhadron=0, pdotqhat=0, Q2i=0;
  Float_t Wi=0, ti=0, phipqi=0, phicm=0, thetapqi=0;
  Float_t saghai=0, factor=0;

  Float_t ppi=0, phi=0, sigcent=0, zi=0, pt2i=0, xbji=0;
  Float_t thqi=0, sighad=0, jacobian=0, centjac=0, xfermi=0;
  Float_t corrsing=0, PmPar=0, PmPer=0, PmOop=0, sigcc=0, radphot=0;
  Float_t M_recoil=0, MM2=0;

  // kinematics for recon
  Float_t had_mom_mag;
  Double_t kf=0, ki=0, Pf=0;
  Double_t Ep=0, En=0, W2=0, fScatAngle=0;

  Double_t th_pq=0; //Polar angle of detected particle with q
  Double_t th_nq=0; //Polar angle of recoil system with q (rad)
  Double_t ph_pq=0; //Azimuth angle of detected particle with q
  Double_t ph_nq=0; //Azimuth of recoil system with scattering plane (rad)

  // rotations / four-vectors
  TRotation fToLabRot; // Rotation matrix from TRANSPORT to lab
  TRotation rot_to_q; // Rotation matrix from +z to +q
  TVector3 Pf_vec, kf_vec, bq, xq, p_miss;

  TLorentzVector fP0; // Beam e- 4-momentum
  TLorentzVector fP1; // Scattered e- 4-momentum
  TLorentzVector fA; // Target 4-momentum
  TLorentzVector fA1; // Final system 4-momentum
  TLorentzVector fQ; // Momentum transfer 4-vector
  TLorentzVector fX; // Detected secondary particle 4-momentum (GeV)
  TLorentzVector fB; // recoil system 4-momentum (GeV)
  TLorentzVector fMp; // Stationary proton 4-mom
  TLorentzVector fMp1; // fMp1=fMp+fQ = Outgoing proton 4-mom

  // rho specific variables
  Float_t Mrho=0, Thrho=0;

  // constants (GeV)
  const Double_t MP  = 0.938272;
  const Double_t MD  = 1.87561;
  const Double_t MN  = 0.939566;
  const Double_t me  = 0.00051099;
  const Double_t mk  = 0.493677;
  const Double_t mpi = 0.139570;
  const Double_t MAL  = 25.131710;
  Double_t tgt_mass  = MP;
  Double_t hadron_mass;

}; // class

// Optional ROOT entry-point wrapper
//void run_recon_hcana(const char* stem, const char* reaction);

#endif
