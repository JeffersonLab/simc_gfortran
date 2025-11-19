// run command: root -l -q -b 'recon_hcana.C+("filename_stem","reaction", "hadron_type", Earm_HMS)'
// example: root -l -b -q 'recon_hcana.C+("heep_12p47deg_7p07gev_hyd_rsidis","heep", "mk", kFALSE)'
// example: root -l -b -q 'recon_hcana.C+("coin_7p87deg_3p632gev_hyd_rsidis","sidis", "mpi")'
// example: root -l -b -q 'recon_hcana.C+("rho_kin37_q24p4_x0p44_z0p67_thpi12p9_pip_lh2","rho")'
// example: root -l -b -q 'recon_hcana.C+("delta_kin37_q24p4_x0p44_z0p67_thpi12p9_pip_lh2","delta")'
// example: root -l -b -q 'recon_hcana.C+("exclusive_kin37_q24p4_x0p44_z0p67_thpi12p9_pip_lh2","exclusive")'

// hadron_type selects final hadrons mass, can be mass of {kaon,pion}. If not passed this argument, mpi will be seleted. For heep, whatever the input, proton mass will be selected which is hard-coded.
// Earm_HMS is either kTRUE or kFALSE. It selects electron arm. If not passed, it will be kTRUE. For heep, electron arm is SHMS, so kFALSE should be passed for heep.
// Don't bother about progress bar not completing 100%. It's computing (i%1000) where i is nentries with nentries = tree->GetEntries();

#include "recon_hcana.h"
using namespace std;

// ---------------- main function ----------------
recon_hcana::recon_hcana(TString filename,
				TString reaction_str,
				TString hadron_type,
				Bool_t Earm_HMS) {
  buildFileName(filename);
  ReadReaction(reaction_str);
  fEarm_HMS = Earm_HMS;

  // Determine hadron mass
  if (hadron_type == "mk") hadron_mass = mk;
  else if (hadron_type == "mpi") hadron_mass = mpi; // Default, if not passed this argument
  else {std::cerr << "Warning: unknown hadron_type '" << hadron_type
                  << "'. Defaulting to pion mass." << std::endl;
    hadron_mass = mpi;
  }

  // parse .hist (plain text)
  simc_nevents    = std::stoi( split( FindString("Ngen",   InSIMCHistname)[0], '=' )[1] );
  simc_normfactor = std::stod( split( FindString("normfac",InSIMCHistname)[0], '=' )[1] );
  Ein             = std::stod( split( FindString("Ebeam",  InSIMCHistname)[0], '=' )[1] ); // MeV

  {
    auto mom_line = split( FindString("momentum", InSIMCHistname)[0], '=' )[1];
    auto ang_line = split( FindString("angle",    InSIMCHistname)[0], '=' )[1];
    auto moms     = num_split(mom_line);
    auto angs     = num_split(ang_line);
    kf0  = moms.at(0);                // MeV
    Pf0  = moms.at(1);                // MeV
    // Depending on electron arm, choose angles
    if(fEarm_HMS){//electron in HMS (-angle), hadron in SHMS (+angle)
      e_th = -angs.at(0);               // deg; HMS sign
      h_th =  angs.at(1);               // deg
    } else{//electron in SHMS (+angle), hadron in SHMS (-angle)
      e_th =  angs.at(0);               // deg
      h_th = -angs.at(1);               // deg; HMS sign
    }
    e_ph = 0.0; h_ph = 0.0;
  }

  cout <<"nevents,normfac,Ein,kf0,Pf0,e_th,h_th: "<<simc_nevents<<","<<simc_normfactor<<","<<Ein<<","<<kf0<<","<<Pf0<<","<<e_th<<","<<h_th<<endl;

  fin  = TFile::Open(InSIMCRootname, "READ");
  if (!fin || fin->IsZombie()) throw std::runtime_error(("Cannot open "+string(InSIMCRootname)).c_str());
  fout = TFile::Open(OutSIMCRootname, "RECREATE");
  if (!fout || fout->IsZombie()) throw std::runtime_error(("Cannot create "+string(OutSIMCRootname)).c_str());

  // Read, Compute, Write
  ReadTree(reaction);
  EventLoop();
  WriteTree();
}

recon_hcana::~recon_hcana() {
  if (fout) { fout->Close(); delete fout; fout=nullptr; }
  if (fin)  { fin->Close();  delete fin;  fin =nullptr; }
}

// --------------- unified reader ----------------
void recon_hcana::ReadTree(const TString& reaction) {
  tree = dynamic_cast<TTree*>(fin->Get("h10"));
  if (!tree) throw std::runtime_error("input tree h10 missing");

  fout->cd();
  newTree = new TTree("h10", "hcana-style recomputed kinematics");
  newTree->SetAutoSave(0);

  nentries = tree->GetEntries();

  // reaction-specific read branches
  if (reaction=="heep") {
    CommonTree();

    tree->SetBranchAddress("corrsing",&corrsing);
    tree->SetBranchAddress("Pmx", &Pmx_simc);
    tree->SetBranchAddress("Pmy", &Pmy_simc);
    tree->SetBranchAddress("Pmz", &Pmz_simc);
    tree->SetBranchAddress("PmPar",&PmPar);
    tree->SetBranchAddress("PmPer",&PmPer);
    tree->SetBranchAddress("PmOop",&PmOop);
    tree->SetBranchAddress("radphot",&radphot);
    tree->SetBranchAddress("sigcc",&sigcc);
  }
    else {
    CommonTree();
    SRDECommonTree(); // SRDE = SidisRhoDeltaExclusive

    // Added if(), so that ROOT doesn't shout warnings if branches are missing.
    if (tree->GetBranch("ppi"))  tree->SetBranchAddress("ppi", &ppi);
    if (tree->GetBranch("sigcent"))  tree->SetBranchAddress("sigcent",   &sigcent);
    if (tree->GetBranch("z"))  tree->SetBranchAddress("z",  &z_simc);
    if (tree->GetBranch("zi"))  tree->SetBranchAddress("zi",  &zi);
    if (tree->GetBranch("pt2"))  tree->SetBranchAddress("pt2",  &pt2_simc);
    if (tree->GetBranch("pt2i"))  tree->SetBranchAddress("pt2i",  &pt2i);
    if (tree->GetBranch("xbj"))  tree->SetBranchAddress("xbj",  &xbj_simc);
    if (tree->GetBranch("xbji"))  tree->SetBranchAddress("xbji",  &xbji);
    if (tree->GetBranch("thqi"))  tree->SetBranchAddress("thqi",  &thqi);
    if (tree->GetBranch("sighad"))  tree->SetBranchAddress("sighad",  &sighad);
    if (tree->GetBranch("jacobian"))  tree->SetBranchAddress("jacobian",  &jacobian);
    if (tree->GetBranch("centjac"))  tree->SetBranchAddress("centjac",  &centjac);
    if (tree->GetBranch("xfermi"))  tree->SetBranchAddress("xfermi",  &xfermi);

    if (tree->GetBranch("mmnuc"))  tree->SetBranchAddress("mmnuc",&mmnuc);
    if (tree->GetBranch("phad"))  tree->SetBranchAddress("phad",&phad);
    if (tree->GetBranch("pmpar"))  tree->SetBranchAddress("pmpar",&pmpar);
    if (tree->GetBranch("pmper"))  tree->SetBranchAddress("pmper",&pmper);
    if (tree->GetBranch("pmoop"))  tree->SetBranchAddress("pmoop",&pmoop);
    if (tree->GetBranch("sigcm"))  tree->SetBranchAddress("sigcm",&sigcm);
    if (tree->GetBranch("pdotqhat"))  tree->SetBranchAddress("pdotqhat",&pdotqhat);
    if (tree->GetBranch("Q2i"))  tree->SetBranchAddress("Q2i",&Q2i);
    if (tree->GetBranch("Wi"))  tree->SetBranchAddress("Wi",&Wi);
    if (tree->GetBranch("ti"))  tree->SetBranchAddress("ti",&ti);

    // rho specific branches
    if (tree->GetBranch("Mrho"))  tree->SetBranchAddress("Mrho",   &Mrho);
    if (tree->GetBranch("Thrho"))  tree->SetBranchAddress("Thrho",   &Thrho);
  }

  // reaction specific write branches
  if (reaction=="heep") {
    CommonNewTree();

    newTree->Branch("corrsing",&corrsing,"corrsing/F");
    newTree->Branch("Pmx",&Pmx_simc,"Pmx/F");
    newTree->Branch("Pmy",&Pmy_simc,"Pmy/F");
    newTree->Branch("Pmz",&Pmz_simc,"Pmz/F");
    newTree->Branch("Pmx_recon",&Pmx_recon,"Pmx_recon/F");
    newTree->Branch("Pmy_recon",&Pmy_recon,"Pmy_recon/F");
    newTree->Branch("Pmz_recon",&Pmz_recon,"Pmz_recon/F");
    newTree->Branch("PmPar",&PmPar,"PmPar/F");
    newTree->Branch("PmPer",&PmPer,"PmPer/F");
    newTree->Branch("PmOop",&PmOop,"PmOop/F");
    newTree->Branch("fry",&fry,"fry/F");
    newTree->Branch("radphot",&radphot,"radphot/F");
    newTree->Branch("sigcc",&sigcc,"sigcc/F");
  }
  else {
    CommonNewTree();
    SRDECommonNewTree(); // SRDE = SidisRhoDeltaExclusive

    newTree->Branch("ppi", &ppi, "ppi/F");
    newTree->Branch("sigcent", &sigcent, "sigcent/F");
    newTree->Branch("z",  &z_simc,  "z/F");
    newTree->Branch("z_recon",  &z_recon,  "z_recon/F");
    newTree->Branch("zi",  &zi,  "zi/F");
    newTree->Branch("pt2",  &pt2_simc,  "pt2/F");
    newTree->Branch("pt2_recon",  &pt2_recon,  "pt2_recon/F");
    newTree->Branch("pt_recon",  &pt_recon,  "pt_recon/F");
    newTree->Branch("ptx_recon",  &ptx_recon,  "ptx_recon/F");
    newTree->Branch("pty_recon",  &pty_recon,  "pty_recon/F");
    newTree->Branch("pt2i",  &pt2i,  "pt2i/F");
    newTree->Branch("xbj",  &xbj_simc,  "xbj/F");
    newTree->Branch("xbj_recon",  &xbj_recon,  "xbj_recon/F");
    newTree->Branch("xbji",  &xbji,  "xbji/F");
    newTree->Branch("thqi",  &thqi,  "thqi/F");
    newTree->Branch("sighad",  &sighad,  "sighad/F");
    newTree->Branch("jacobian",  &jacobian,  "jacobian/F");
    newTree->Branch("centjac",  &centjac,  "centjac/F");
    newTree->Branch("xfermi",  &xfermi,  "xfermi/F");

    newTree->Branch("mmnuc",   &mmnuc,   "mmnuc/F");
    newTree->Branch("phad",   &phad,   "phad/F");
    newTree->Branch("pmpar",   &pmpar,   "pmpar/F");
    newTree->Branch("pmper",   &pmper,   "pmper/F");
    newTree->Branch("pmoop",   &pmoop,   "pmoop/F");
    newTree->Branch("sigcm",   &sigcm,   "sigcm/F");
    newTree->Branch("pdotqhat",   &pdotqhat,   "pdotqhat/F");
    newTree->Branch("Q2i",   &Q2i,   "Q2i/F");
    newTree->Branch("Wi",   &Wi,   "Wi/F");
    newTree->Branch("ti",   &ti,   "ti/F");
    // rho specific branches
    newTree->Branch("Mrho",   &Mrho,   "Mrho/F");
    newTree->Branch("Thrho",   &Thrho,   "Thrho/F");
  }
}

void recon_hcana::CommonTree(){
    tree->SetBranchAddress("hsdelta", &hsdelta);
    tree->SetBranchAddress("hsyptar", &hsyptar);
    tree->SetBranchAddress("hsxptar", &hsxptar);
    tree->SetBranchAddress("hsytar",  &hsytar);
    tree->SetBranchAddress("hsxfp",  &hsxfp);
    tree->SetBranchAddress("hsxpfp",  &hsxpfp);
    tree->SetBranchAddress("hsyfp",  &hsyfp);
    tree->SetBranchAddress("hsypfp",  &hsypfp);
    tree->SetBranchAddress("hsdeltai",  &hsdeltai);
    tree->SetBranchAddress("hsxptari",  &hsxptari);
    tree->SetBranchAddress("hsyptari",  &hsyptari);
    tree->SetBranchAddress("hsytari",  &hsytari);

    tree->SetBranchAddress("ssdelta", &ssdelta);
    tree->SetBranchAddress("ssyptar", &ssyptar);
    tree->SetBranchAddress("ssxptar", &ssxptar);
    tree->SetBranchAddress("ssytar",  &ssytar);
    tree->SetBranchAddress("ssxfp",  &ssxfp);
    tree->SetBranchAddress("ssxpfp",  &ssxpfp);
    tree->SetBranchAddress("ssyfp",  &ssyfp);
    tree->SetBranchAddress("ssypfp",  &ssypfp);
    tree->SetBranchAddress("ssdeltai",  &ssdeltai);
    tree->SetBranchAddress("ssxptari",  &ssxptari);
    tree->SetBranchAddress("ssyptari",  &ssyptari);
    tree->SetBranchAddress("ssytari",  &ssytari);

    tree->SetBranchAddress("Weight",  &Weight);
    tree->SetBranchAddress("q",  &q_simc);
    tree->SetBranchAddress("nu", &nu_simc);
    tree->SetBranchAddress("Q2", &Q2_simc);
    tree->SetBranchAddress("W",  &W_simc);
    tree->SetBranchAddress("epsilon", &epsilon_simc);
    tree->SetBranchAddress("Em",      &Em_simc);
    tree->SetBranchAddress("Pm",      &Pm_simc);
    tree->SetBranchAddress("thetapq", &thetapq_simc);
    tree->SetBranchAddress("phipq",   &phipq_simc);
  }

void recon_hcana::CommonNewTree(){
    newTree->Branch("Ngen",    &simc_nevents,    "Ngen/I");
    newTree->Branch("normfac", &simc_normfactor, "normfac/D");
    newTree->Branch("Weight",   &Weight,   "Weight/F");
    newTree->Branch("fWeight", &fWeight,         "fWeight/F");

    newTree->Branch("hsdelta", &hsdelta, "hsdelta/F");
    newTree->Branch("hsyptar", &hsyptar, "hsyptar/F");
    newTree->Branch("hsxptar", &hsxptar, "hsxptar/F");
    newTree->Branch("hsytar",  &hsytar,  "hsytar/F");
    newTree->Branch("hsxfp",  &hsxfp,  "hsxfp/F");
    newTree->Branch("hsxpfp",  &hsxpfp,  "hsxpfp/F");
    newTree->Branch("hsyfp",  &hsyfp,  "hsyfp/F");
    newTree->Branch("hsypfp",  &hsypfp,  "hsypfp/F");
    newTree->Branch("hsdeltai",  &hsdeltai,  "hsdeltai/F");
    newTree->Branch("hsxptari",  &hsxptari,  "hsxptari/F");
    newTree->Branch("hsyptari",  &hsyptari,  "hsyptari/F");
    newTree->Branch("hsytari",  &hsytari,  "hsytari/F");

    newTree->Branch("ssdelta", &ssdelta, "ssdelta/F");
    newTree->Branch("ssyptar", &ssyptar, "ssyptar/F");
    newTree->Branch("ssxptar", &ssxptar, "ssxptar/F");
    newTree->Branch("ssytar",  &ssytar,  "ssytar/F");
    newTree->Branch("ssxfp",  &ssxfp,  "ssxfp/F");
    newTree->Branch("ssxpfp",  &ssxpfp,  "ssxpfp/F");
    newTree->Branch("ssyfp",  &ssyfp,  "ssyfp/F");
    newTree->Branch("ssypfp",  &ssypfp,  "ssypfp/F");
    newTree->Branch("ssdeltai",  &ssdeltai,  "ssdeltai/F");
    newTree->Branch("ssxptari",  &ssxptari,  "ssxptari/F");
    newTree->Branch("ssyptari",  &ssyptari,  "ssyptari/F");
    newTree->Branch("ssytari",  &ssytari,  "ssytari/F");

    newTree->Branch("q",  &q_simc,        "q/F");
    newTree->Branch("nu", &nu_simc,       "nu/F");
    newTree->Branch("Q2", &Q2_simc,  "Q2/F");
    newTree->Branch("W",  &W_simc,   "W/F");
    newTree->Branch("epsilon",  &epsilon_simc,  "epsilon/F");
    newTree->Branch("Em",       &Em_simc,       "Em/F");
    newTree->Branch("Pm",       &Pm_simc,       "Pm/F");
    newTree->Branch("thetapq",  &thetapq_simc,  "thetapq/F");
    newTree->Branch("phipq",    &phipq_simc,    "phipq/F");

    newTree->Branch("q_recon", &q_recon, "q_recon/F");
    newTree->Branch("nu_recon", &nu_recon, "nu_recon/F");
    newTree->Branch("Q2_recon", &Q2_recon, "Q2_recon/F");
    newTree->Branch("W_recon",  &W_recon,  "W_recon/F");
    newTree->Branch("epsilon_recon",  &epsilon_recon,  "epsilon_recon/F");
    newTree->Branch("Em_recon", &Em_recon, "Em_recon/F");
    newTree->Branch("Pm_recon", &Pm_recon, "Pm_recon/F");
    newTree->Branch("thetapq_recon", &thetapq_recon, "thetapq_recon/F");
    newTree->Branch("phipq_recon", &phipq_recon, "phipq_recon/F");
  }

void recon_hcana::SRDECommonTree(){ // SRDE = SidisRhoDeltaExclusive
    tree->SetBranchAddress("missmass", &missmass_simc);
    tree->SetBranchAddress("t",        &t);
    tree->SetBranchAddress("fry",      &fry);
    tree->SetBranchAddress("radphot",  &radphot);
    tree->SetBranchAddress("siglab",   &siglab);
    tree->SetBranchAddress("decdist",  &decdist);
    tree->SetBranchAddress("Mhadron",  &Mhadron);
    tree->SetBranchAddress("pfermi",   &pfermi);
    tree->SetBranchAddress("phipqi",  &phipqi);
  }

void recon_hcana::SRDECommonNewTree(){ // SRDE = SidisRhoDeltaExclusive
    newTree->Branch("missmass", &missmass_simc, "missmass/F");
    newTree->Branch("missmass_recon", &missmass_recon, "missmass_recon/F");
    newTree->Branch("t",        &t,        "t/F");
    newTree->Branch("fry",      &fry,      "fry/F");
    newTree->Branch("radphot",  &radphot,  "radphot/F");
    newTree->Branch("siglab",   &siglab,   "siglab/F");
    newTree->Branch("decdist",  &decdist,  "decdist/F");
    newTree->Branch("Mhadron",  &Mhadron,  "Mhadron/F");
    newTree->Branch("pfermi",   &pfermi,   "pfermi/F");
    newTree->Branch("phipqi",  &phipqi,  "phipqi/F");
  }


// --------------- event loop ----------------
void recon_hcana::EventLoop() {
  // MeV -> GeV for recon
  Ein /= 1000.0; //Ebeam from simc .hist
  kf0 /= 1000.0; //E-arm central momentum
  Pf0 /= 1000.0; //P-arm central momentum

  for (Int_t i=0;i<nentries;i++) {
    if (reaction != "heep"){
    //Set ErrorVal to those branches which are missing
    ppi = sigcent = z_simc = z_recon = ErrorVal;
    zi = pt2_simc = pt2_recon = pt_recon = ErrorVal;
    ptx_recon = pty_recon = pt2i = xbj_simc = ErrorVal;
    xbji = thqi = sighad = jacobian = ErrorVal;
    centjac = xfermi = xbj_recon = ErrorVal;

    mmnuc = phad = pmpar = pmper = ErrorVal;
    pmoop = sigcm = pdotqhat = Q2i = ErrorVal;
    Wi = ti = ErrorVal;

    Mrho = Thrho = ErrorVal;
    }

    // Progress bar
    if(i%1000==0) {
      int barWidth = 25;
      progress = ((float)i/(float)nentries);
      cout << "[";
      float pos = barWidth * progress;
      for (float i = 0.; i < barWidth; ++i) {
	if (i < pos) cout << "=";
	else if (i == pos) cout << ">";
	else cout << " ";
      }
      cout << "] " << int(progress * 100.0) << " %\r";
      cout.flush();
    }

    tree->GetEntry(i);

    if (fEarm_HMS) {// electron in HMS, hadron in SHMS (default behavior)
      edelta = hsdelta; exptar = hsxptar; eyptar = hsyptar;
      hdelta = ssdelta; hxptar = ssxptar; hyptar = ssyptar;
    } else {// electron in SHMS, hadron in HMS (heep case)
      edelta = ssdelta; exptar = ssxptar; eyptar = ssyptar;
      hdelta = hsdelta; hxptar = hsxptar; hyptar = hsyptar;
    }

    // momentum from delta: delta = 100*(mom_track - mom_central)/mom_central
    kf = kf0 * (1.0 + edelta/100.0); //E-arm mom_track
    Pf = Pf0 * (1.0 + hdelta/100.0); //P-arm mom_track

    // electron 4-vectors
    ki = std::sqrt(Ein*Ein - me*me); //initial e- mom mag
    SetCentralAngles(e_th, e_ph);
    TransportToLab(kf, exptar, eyptar, kf_vec); //calculates: TVector3 kf_vec
    fP0.SetXYZM(0,0,ki,me); //set incoming e- 4-mom
    fP1.SetXYZM(kf_vec.X(), kf_vec.Y(), kf_vec.Z(), me); //set outgoing e- 4-mom in lab frame
    fA.SetXYZM(0,0,0,tgt_mass); //set stationary target 4-mom

    fQ = fP0 - fP1; //set (small q) q = incoming - outgoing
    fScatAngle = fP0.Angle(fP1.Vect()); //calculates: scattering angle between incoming and outgoing e-

    // recon
    fMp.SetXYZM(0,0,0,MP); //set stationary proton 4-mom
    fMp1 = fMp + fQ; //set outgoing proton 4-mom
    W2 = fMp1.M2(); //calculates: invariant mass-squared = E^2 - P^2 for outgoing proton
    q_recon = fQ.P(); //P() calculates 3-mom mag of fQ.
    Q2_recon = -fQ.M2(); //Q^2 = negative of transferred 4-mom squared
    W_recon  = (W2>0)? std::sqrt(W2):0.0; //invariant mass
    epsilon_recon = 1.0 / ( 1.0 + 2.0*q_recon*q_recon/Q2_recon*TMath::Power( TMath::Tan(fScatAngle/2.0), 2.0 )); //photon polarization vectors magnitude. See Halzen and Martin eqn 8.57

    // detected/outgoing hadron 4-vector
    SetCentralAngles(h_th, h_ph);
    TransportToLab(Pf, hxptar, hyptar, Pf_vec); //calculates: TVector3 Pf_vec. Outgoing hadron.
    //if heep, set outgoing hadron mass proton. otherwise pion (for rsidis)
    double mass_h = (reaction=="heep" ? double(MP) : hadron_mass);
    fX.SetVectM(Pf_vec, mass_h); //set outgoing hadron 4-mom

    fA1 = fA + fQ; //final target 4-mom
    fB  = fA1 - fX; //recoil 4-mom

    nu_recon = fQ.E(); //nu = E - E-prime. .E()->energy component of fQ
    Em_recon = nu_recon + fA.M() - fX.E(); //.M()->invariant mass. Em=recoil energy. Missing energy.

    //--------Rotate the recoil system from +z to +q-------
    // Angles of X and B wrt q-vector
    // xq and bq are the 3-momentum vectors of X and B expressed in
    // the coordinate system where q is the z-axis and the x-axis
    // lies in the scattering plane (defined by q and e') and points
    // in the direction of e', so the out-of-plane angle lies within
    // -90<phi_xq<90deg if X is detected on the downstream/forward side of q
    rot_to_q.SetZAxis( fQ.Vect(), fP1.Vect()).Invert();
    //.SetZAxis() sets z-axis along 1st_argument=fQ and normalize it to unit vector, then internally
    //calculates the cross product of new z-axis and 2nd argument (an arbitrary vector) fP1. The
    //resulting vector will be used to define new y-axis (typically). A second cross product is then
    //calculated to get x-axis. All axis vectors get normalized to unity.
    //{fQ,fP1}.Vect() are lab frame vector according to their declaration. They are being used to
    //create a rotation matrix, say R, to go from q to lab frame. So,
    //R = rot_to_q.SetZAxis( fQ.Vect(), fP1.Vect()). This (R) matrix elements are q-frames unit vectors
    //in lab frame. So if you do R*any_vector_in_q_frame, you'll get that vector in lab-frame.
    //Therefore, we take the .Invert() to get transformation for lab to q-frame.

    xq = fX.Vect(); //3-mom of fX=outgoing hadron 4-mom
    bq = fB.Vect(); //3-mom of fB=reoil 4-mom

    xq *= rot_to_q; //matrix-vetor multiplication. New_xq = rot_to_q * Old_xq
    bq *= rot_to_q; //similar for recoil.
    //so new xq, bq are in q-frame

    //Calculate Angles of q relative to x(detected hadron) and b(recoil neutron)
    th_pq = xq.Theta(); //spherical polar angle of 3D vector xq.
    ph_pq = xq.Phi(); //spherical azimuthal (phi) angle of 3D vector xq
    th_nq = bq.Theta(); //similarly for recoil (neutron) polar
    ph_nq = bq.Phi(); //similarly for recoil (neutron) azimuthal
    //th_pq (th_nq) defines angle between hadron(recoil) w.r.t q bcz q is z axis.

    if(reaction == "sidis" || reaction == "rho"){
    //calculate z,xbj,pt,ptx,pty,pt2
    had_mom_mag = Pf_vec.Mag();
    z_recon = had_mom_mag/nu_recon;
    xbj_recon = Q2_recon/(2*MP*nu_recon);
    pt_recon = TMath::Sqrt(TMath::Power(had_mom_mag,2)*(1.0-TMath::Power(TMath::Cos(th_pq),2)));
    ptx_recon = pt_recon*TMath::Cos(ph_pq);
    pty_recon = pt_recon*TMath::Sin(ph_pq);
    pt2_recon = TMath::Power(pt_recon,2);
    }

    // Explicit naming of angle variables
    thetapq_recon = th_pq; //angle between hadron and q which is the spherical polar angle
    phipq_recon = ph_pq; //angle of the projection of xq (hadron) on x-y plane w.r.t x-axis

    p_miss = -bq; //negative of TVector3 bq (recoil).This is missing 3-mom. Bcz 3-mom is conserved
    //in heep.So, P_init=P_final =>P_init=xq + negative of recoil. Negative of recoil is missing 3-mom.

    //Missing Momentum Components in the q-frame
    Pm_recon = p_miss.Mag();

    // Redefine variables
    Pmx_recon = p_miss.X(); //in-plane perpendicular component to +z
    Pmy_recon = p_miss.Y(); //out-of-plane component (Oop)
    Pmz_recon = p_miss.Z(); //parallel component to +z

    M_recoil = fB.M(); //recoil mass (missing mass)

    missmass_recon = sqrt(abs((Em_recon*Em_recon)-(Pm_recon*Pm_recon)));

    MM2 = missmass_recon * missmass_recon;

    //s = (fQ+fA).M2();
    //t = (fQ-fX).M2();
    //u = (fQ-fB).M2();

    // reconstructed Weight branch
    fWeight = (simc_nevents>0) ? float_t(Weight * simc_normfactor / simc_nevents) : 0.0f;

    newTree->Fill();
  }
}

// --------------- write ----------------
void recon_hcana::WriteTree() {
  fout->cd();
  newTree->Write("h10");
}

// --------------- aux  ----------------
//convert geographical to spherical angle. Units are rad.
void recon_hcana::GeoToSph(Double_t th_geo, Double_t ph_geo, Double_t& th_sph, Double_t& ph_sph) {
  static const Double_t twopi = 2.0*TMath::Pi();
  Double_t ct = cos(th_geo), cp = cos(ph_geo);
  Double_t tmp = ct*cp;
  th_sph = acos(tmp);
  tmp = sqrt(1.0 - tmp*tmp);
  ph_sph = (fabs(tmp) < 1e-6 ) ? 0.0 : acos( sqrt(1.0-ct*ct)*cp/tmp );
  if( th_geo/twopi-floor(th_geo/twopi) > 0.5 ) ph_sph = TMath::Pi() - ph_sph;
  if( ph_geo/twopi-floor(ph_geo/twopi) > 0.5 ) ph_sph = -ph_sph;
}

void recon_hcana::SetCentralAngles(Double_t th_cent=0, Double_t ph_cent=0) {
  Double_t fThetaGeo = TMath::DegToRad()*th_cent;
  Double_t fPhiGeo   = TMath::DegToRad()*ph_cent;
  Double_t fThetaSph, fPhiSph, fSinThGeo, fCosThGeo, fSinPhGeo, fCosPhGeo;
  GeoToSph(fThetaGeo, fPhiGeo, fThetaSph, fPhiSph);
  fSinThGeo = TMath::Sin( fThetaGeo ); fCosThGeo = TMath::Cos( fThetaGeo );
  fSinPhGeo = TMath::Sin( fPhiGeo );   fCosPhGeo = TMath::Cos( fPhiGeo );
  Double_t st = TMath::Sin( fThetaSph ), ct = TMath::Cos( fThetaSph );
  Double_t sp = TMath::Sin( fPhiSph ),   cp = TMath::Cos( fPhiSph );
  Double_t norm = TMath::Sqrt(ct*ct + st*st*cp*cp);
  TVector3 nx( st*st*sp*cp/norm, -norm,          st*ct*sp/norm );
  TVector3 ny( ct/norm,          0.0,            -st*cp/norm   );
  TVector3 nz( st*cp,            st*sp,          ct            );
  fToLabRot.SetToIdentity().RotateAxes(nx,ny,nz);
}

void recon_hcana::TransportToLab(Double_t p, Double_t xptar, Double_t yptar, TVector3& pvect) {
  TVector3 v(xptar, yptar, 1.0);
  v *= p / TMath::Sqrt(1.0 + xptar*xptar + yptar*yptar);
  pvect = fToLabRot * v;
}

// --- text utils ---
std::vector<std::string> recon_hcana::FindString(TString keyword, TString fname) {
  ifstream ifile(fname.Data());
  std::vector<std::string> lines; std::string line;
  while (std::getline(ifile, line))
    if (line.find(keyword.Data()) != std::string::npos) lines.push_back(line);
  return lines;
}
std::vector<std::string> recon_hcana::split(std::string str, char del){
  int p = -1; for (int i=0;i<(int)str.size();++i) if (str[i]==del) p=i;
  std::string L,R; for (int i=0;i<(int)str.size();++i){ if(i<p) L.push_back(str[i]); else if(i>p) R.push_back(str[i]); }
  return {L,R};
}
std::vector<float> recon_hcana::num_split(std::string s){ std::istringstream ss(s); std::vector<float> v; float x; while (ss>>x) v.push_back(x); return v; }
std::string recon_hcana::getString(char x){ return std::string(1,x); }

// ROOT wrapper
//void run_recon_hcana(const char* stem, const char* reaction){ recon_hcana r(stem, reaction); }
