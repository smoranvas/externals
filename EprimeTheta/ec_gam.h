//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug  9 11:57:22 2017 by ROOT version 6.10/02
// from TTree ec_gam/TreeGSIM
// found on file: carbon_id1.root
//////////////////////////////////////////////////////////


#ifndef ec_gam_h
#define ec_gam_h


#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TVirtualIndex.h>
#include <TCanvas.h>
#include <iostream>             // std::cout, std::endl
#include <fstream>              // std::ifstream
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TH2F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TPostScript.h>
#include <TPaveText.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <TStyle.h>
#include <TLatex.h>

using namespace std;
const float qmin= 1.;
const float qmax= 4.1;
const float vmin= 2.2;
const float vmax= 4.25;

const int nqbins=3;
const int nvbins=3;

const float q_values[4]= {1.0, 1.33, 1.76, 4.1 };
const float v_values[4]= {2.2, 3.2, 3.73, 4.25 };
const float center_q[3] = { 1.165 , 1.55, 2.93 };
const float center_v[3] = { 2.7, 3.46, 3.99 };

const float Ebeam=5.014;

class ec_gam : public TSelector {
public :
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   TFile *NewFile;
   TTree *TreeGam;

   Int_t AvailableEvents, Nevta, Nevt_found;
   TString TargetChoice;

   Float_t Eprime, Theta;
   Float_t mean_theta[nqbins][nvbins];
   Float_t mean_Eprime[nqbins][nvbins];

   Int_t qbin, vbin;
   TH2F *h_etheta;
   TH2F *h_etheta_b[nqbins][nvbins];
   TF1 *Q2lines; 
   TLine *nu;

   Float_t         trig_ql;
   Float_t         tr_time;
   UInt_t          evntid;
   Int_t           trig_sect;
   Int_t           sect;
   Float_t         trig_mom;
   Float_t         trig_theta;
   Float_t         trig_phi;
   Float_t         trig_phisect;
   Float_t         trig_nphe;
   Float_t         trig_ecin;
   Float_t         trig_ecout;
   Float_t         trig_ectot;
   Float_t         trig_ece;
   Float_t         trig_vz;
   Float_t         trig_vx;
   Float_t         trig_vy;
   Float_t         trig_vx_corr;
   Float_t         trig_vy_corr;
   Float_t         trig_vz_corr;
   Float_t         trig_px;
   Float_t         trig_py;
   Float_t         trig_pz;
   Float_t         trig_ecx;
   Float_t         trig_ecy;
   Float_t         trig_ecz;
   Float_t         trig_ecu;
   Float_t         trig_ecv;
   Float_t         trig_ecw;
   Float_t         trig_ect;
   Float_t         trig_ecl;
   Float_t         trig_sct;
   Float_t         trig_scl;
   Int_t           trig_sc_stat;
   Float_t         trig_cctheta;
   Float_t         trig_ccphi;
   Float_t         trig_tl1x;
   Float_t         trig_tl1y;
   Float_t         trig_tl1z;
   Float_t         trig_tl1cx;
   Float_t         trig_tl1cy;
   Float_t         trig_tl1cz;
   Float_t         trig_phil1;
   Float_t         trig_thetal1;
   Float_t         trig_W;
   Float_t         trig_Q2;
   Float_t         trig_nu;
   Float_t         trig_y;
   Int_t           target;
   Float_t         trig_transv_pos;
   Int_t           LorenzoFlag;
   Bool_t          NpheCutSect;
   Int_t           trig_dcfid_flag;
   Int_t           MirrorCode;
   Int_t           nphot;
   Float_t         gam_Eec[8];   //[nphot]
   Float_t         gam_bet[8];   //[nphot]
   Float_t         gam_Pxtrue[8];   //[nphot]
   Float_t         gam_Pytrue[8];   //[nphot]
   Float_t         gam_Pztrue[8];   //[nphot]
   Float_t         gam_E[8];   //[nphot]
   Float_t         gam_Px[8];   //[nphot]
   Float_t         gam_Py[8];   //[nphot]
   Float_t         gam_Pz[8];   //[nphot]
   Float_t         gam_Etrue[8];   //[nphot]
   Float_t         gam_Etrue_Corr[8];   //[nphot]
   Float_t         gam_Pxtrue_Corr[8];   //[nphot]
   Float_t         gam_Pytrue_Corr[8];   //[nphot]
   Float_t         gam_Pztrue_Corr[8];   //[nphot]
   Float_t         gam_ecx[8];   //[nphot]
   Float_t         gam_ecy[8];   //[nphot]
   Float_t         gam_ecz[8];   //[nphot]
   Float_t         gam_ecu[8];   //[nphot]
   Float_t         gam_ecv[8];   //[nphot]
   Float_t         gam_ecw[8];   //[nphot]
   Float_t         gam_ect[8];   //[nphot]
   Float_t         gam_ecpath[8];   //[nphot]
   Float_t         gam_phi[8];   //[nphot]
   Float_t         gam_theta[8];   //[nphot]
   Int_t           npi0;
   Float_t         Pi0z_EVNT[28];   //[npi0]
   Float_t         Pi0mass_EVNT[28];   //[npi0]
   Float_t         opening_angle_evnt[28];   //[npi0]
   Float_t         Pi0energy_EVNT[28];   //[npi0]
   Float_t         Pi0z_ECPB[28];   //[npi0]
   Float_t         Pi0mass_ECPB[28];   //[npi0]
   Float_t         opening_angle_ecpb[28];   //[npi0]
   Float_t         Pi0energy_ECPB[28];   //[npi0]
   Float_t         Pi0mass_Corr[28];   //[npi0]
   Float_t         Pi0z_Corr[28];   //[npi0]
   Float_t         Pi0energy_Corr[28];   //[npi0]
   Float_t         Pi0pheta_Corr[28];   //[npi0]
   Float_t         Pi0phi_Corr[28];   //[npi0]
   Float_t         Pi0angle_gs[28];   //[npi0]
   Float_t         pT2_val_Corr[28];   //[npi0]
   Float_t         phi_val_Corr[28];   //[npi0]
   Float_t         pT1_val[28];   //[npi0]
   Float_t         opening_angle_Corr[28];   //[npi0]
   Int_t           G1[28];   //[npi0]
   Int_t           G2[28];   //[npi0]
   Float_t         eG1Angle[28];   //[npi0]
   Float_t         eG2Angle[28];   //[npi0]
   // List of branches
   TBranch        *b_trig_ql;   //!
   TBranch        *b_tr_time;   //!
   TBranch        *b_evntid;   //!
   TBranch        *b_trig_sect;   //!
   TBranch        *b_sect;   //!
   TBranch        *b_trig_mom;   //!
   TBranch        *b_trig_theta;   //!
   TBranch        *b_trig_phi;   //!
   TBranch        *b_trig_phisect;   //!
   TBranch        *b_trig_nphe;   //!
   TBranch        *b_trig_ecin;   //!
   TBranch        *b_trig_ecout;   //!
   TBranch        *b_trig_ectot;   //!
   TBranch        *b_trig_ece;   //!
   TBranch        *b_trig_vz;   //!
   TBranch        *b_trig_vx;   //!
   TBranch        *b_trig_vy;   //!
   TBranch        *b_trig_vx_corr;   //!
   TBranch        *b_trig_vy_corr;   //!
   TBranch        *b_trig_vz_corr;   //!
   TBranch        *b_trig_px;   //!
   TBranch        *b_trig_py;   //!
   TBranch        *b_trig_pz;   //!
   TBranch        *b_trig_ecx;   //!
   TBranch        *b_trig_ecy;   //!
   TBranch        *b_trig_ecz;   //!
   TBranch        *b_trig_ecu;   //!
   TBranch        *b_trig_ecv;   //!
   TBranch        *b_trig_ecw;   //!
   TBranch        *b_trig_ect;   //!
   TBranch        *b_trig_ecl;   //!
   TBranch        *b_trig_sct;   //!
   TBranch        *b_trig_scl;   //!
   TBranch        *b_trig_sc_stat;   //!
   TBranch        *b_trig_cctheta;   //!
   TBranch        *b_trig_ccphi;   //!
   TBranch        *b_trig_tl1x;   //!
   TBranch        *b_trig_tl1y;   //!
   TBranch        *b_trig_tl1z;   //!
   TBranch        *b_trig_tl1cx;   //!
   TBranch        *b_trig_tl1cy;   //!
   TBranch        *b_trig_tl1cz;   //!
   TBranch        *b_trig_phil1;   //!
   TBranch        *b_trig_thetal1;   //!
   TBranch        *b_trig_W;   //!
   TBranch        *b_trig_Q2;   //!
   TBranch        *b_trig_nu;   //!
   TBranch        *b_trig_y;   //!
   TBranch        *b_target;   //!
   TBranch        *b_trig_transv_pos;   //!
   TBranch        *b_LorenzoFlag;   //!
   TBranch        *b_NpheCutSect;   //!
   TBranch        *b_trig_dcfid_flag;   //!
   TBranch        *b_MirrorCode;   //!
   TBranch        *b_nphot;   //!
   TBranch        *b_gam_Eec;   //!
   TBranch        *b_gam_bet;   //!
   TBranch        *b_gam_Pxtrue;   //!
   TBranch        *b_gam_Pytrue;   //!
   TBranch        *b_gam_Pztrue;   //!
   TBranch        *b_gam_E;   //!
   TBranch        *b_gam_Px;   //!
   TBranch        *b_gam_Py;   //!
   TBranch        *b_gam_Pz;   //!
   TBranch        *b_gam_Etrue;   //!
   TBranch        *b_gam_Etrue_Corr;   //!
   TBranch        *b_gam_Pxtrue_Corr;   //!
   TBranch        *b_gam_Pytrue_Corr;   //!
   TBranch        *b_gam_Pztrue_Corr;   //!
   TBranch        *b_gam_ecx;   //!
   TBranch        *b_gam_ecy;   //!
   TBranch        *b_gam_ecz;   //!
   TBranch        *b_gam_ecu;   //!
   TBranch        *b_gam_ecv;   //!
   TBranch        *b_gam_ecw;   //!
   TBranch        *b_gam_ect;   //!
   TBranch        *b_gam_ecpath;   //!
   TBranch        *b_gam_phi;   //!
   TBranch        *b_gam_theta;   //!
   TBranch        *b_npi0;   //!
   TBranch        *b_Pi0z_EVNT;   //!
   TBranch        *b_Pi0mass_EVNT;   //!
   TBranch        *b_opening_angle_evnt;   //!
   TBranch        *b_Pi0energy_EVNT;   //!
   TBranch        *b_Pi0z_ECPB;   //!
   TBranch        *b_Pi0mass_ECPB;   //!
   TBranch        *b_opening_angle_ecpb;   //!
   TBranch        *b_Pi0energy_ECPB;   //!
   TBranch        *b_Pi0mass_Corr;   //!
   TBranch        *b_Pi0z_Corr;   //!
   TBranch        *b_Pi0energy_Corr;   //!
   TBranch        *b_Pi0pheta_Corr;   //!
   TBranch        *b_Pi0phi_Corr;   //!
   TBranch        *b_Pi0angle_gs;   //!
   TBranch        *b_pT2_val_Corr;   //!
   TBranch        *b_phi_val_Corr;   //!
   TBranch        *b_pT1_val;   //!
   TBranch        *b_opening_angle_Corr;   //!
   TBranch        *b_G1;   //!
   TBranch        *b_G2;   //!
   TBranch        *b_eG1Angle;   //!
   TBranch        *b_eG2Angle;   //!
 
    ec_gam (TTree * tree =0, int Ntoproc=0, TString MyTargetChoice="") {AvailableEvents=Ntoproc; TargetChoice=MyTargetChoice;}
   virtual ~ec_gam() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   void Calculate(void);
   void BookHists(void);
   void WriteHists(void);
   void WriteTxt(void);
   void FillTables(void);
   void MakeBranch(void);

   int Get_q_bin(float q);
   int Get_v_bin(float v);

   ClassDef(ec_gam,0);

};

#endif

#ifdef ec_gam_cxx
void ec_gam::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).
   if (!tree) return;
   fChain = tree;
   //fCurrent = -1;
   fChain->SetMakeClass(1);
   fChain->SetBranchAddress("trig_ql", &trig_ql, &b_trig_ql);
   fChain->SetBranchAddress("tr_time", &tr_time, &b_tr_time);
   fChain->SetBranchAddress("evntid", &evntid, &b_evntid);
   fChain->SetBranchAddress("trig_sect", &trig_sect, &b_trig_sect);
   fChain->SetBranchAddress("sect", &sect, &b_sect);
   fChain->SetBranchAddress("trig_mom", &trig_mom, &b_trig_mom);
   fChain->SetBranchAddress("trig_theta", &trig_theta, &b_trig_theta);
   fChain->SetBranchAddress("trig_phi", &trig_phi, &b_trig_phi);
   fChain->SetBranchAddress("trig_phisect", &trig_phisect, &b_trig_phisect);
   fChain->SetBranchAddress("trig_nphe", &trig_nphe, &b_trig_nphe);
   fChain->SetBranchAddress("trig_ecin", &trig_ecin, &b_trig_ecin);
   fChain->SetBranchAddress("trig_ecout", &trig_ecout, &b_trig_ecout);
   fChain->SetBranchAddress("trig_ectot", &trig_ectot, &b_trig_ectot);
   fChain->SetBranchAddress("trig_ece", &trig_ece, &b_trig_ece);
   fChain->SetBranchAddress("trig_vz", &trig_vz, &b_trig_vz);
   fChain->SetBranchAddress("trig_vx", &trig_vx, &b_trig_vx);
   fChain->SetBranchAddress("trig_vy", &trig_vy, &b_trig_vy);
   fChain->SetBranchAddress("trig_vx_corr", &trig_vx_corr, &b_trig_vx_corr);
   fChain->SetBranchAddress("trig_vy_corr", &trig_vy_corr, &b_trig_vy_corr);
   fChain->SetBranchAddress("trig_vz_corr", &trig_vz_corr, &b_trig_vz_corr);
   fChain->SetBranchAddress("trig_px", &trig_px, &b_trig_px);
   fChain->SetBranchAddress("trig_py", &trig_py, &b_trig_py);
   fChain->SetBranchAddress("trig_pz", &trig_pz, &b_trig_pz);
   fChain->SetBranchAddress("trig_ecx", &trig_ecx, &b_trig_ecx);
   fChain->SetBranchAddress("trig_ecy", &trig_ecy, &b_trig_ecy);
   fChain->SetBranchAddress("trig_ecz", &trig_ecz, &b_trig_ecz);
   fChain->SetBranchAddress("trig_ecu", &trig_ecu, &b_trig_ecu);
   fChain->SetBranchAddress("trig_ecv", &trig_ecv, &b_trig_ecv);
   fChain->SetBranchAddress("trig_ecw", &trig_ecw, &b_trig_ecw);
   fChain->SetBranchAddress("trig_ect", &trig_ect, &b_trig_ect);
   fChain->SetBranchAddress("trig_ecl", &trig_ecl, &b_trig_ecl);
   fChain->SetBranchAddress("trig_sct", &trig_sct, &b_trig_sct);
   fChain->SetBranchAddress("trig_scl", &trig_scl, &b_trig_scl);
   fChain->SetBranchAddress("trig_sc_stat", &trig_sc_stat, &b_trig_sc_stat);
   fChain->SetBranchAddress("trig_cctheta", &trig_cctheta, &b_trig_cctheta);
   fChain->SetBranchAddress("trig_ccphi", &trig_ccphi, &b_trig_ccphi);
   fChain->SetBranchAddress("trig_tl1x", &trig_tl1x, &b_trig_tl1x);
   fChain->SetBranchAddress("trig_tl1y", &trig_tl1y, &b_trig_tl1y);
   fChain->SetBranchAddress("trig_tl1z", &trig_tl1z, &b_trig_tl1z);
   fChain->SetBranchAddress("trig_tl1cx", &trig_tl1cx, &b_trig_tl1cx);
   fChain->SetBranchAddress("trig_tl1cy", &trig_tl1cy, &b_trig_tl1cy);
   fChain->SetBranchAddress("trig_tl1cz", &trig_tl1cz, &b_trig_tl1cz);
   fChain->SetBranchAddress("trig_phil1", &trig_phil1, &b_trig_phil1);
   fChain->SetBranchAddress("trig_thetal1", &trig_thetal1, &b_trig_thetal1);
   fChain->SetBranchAddress("trig_W", &trig_W, &b_trig_W);
   fChain->SetBranchAddress("trig_Q2", &trig_Q2, &b_trig_Q2);
   fChain->SetBranchAddress("trig_nu", &trig_nu, &b_trig_nu);
   fChain->SetBranchAddress("trig_y", &trig_y, &b_trig_y);
   fChain->SetBranchAddress("target", &target, &b_target);
   fChain->SetBranchAddress("trig_transv_pos", &trig_transv_pos, &b_trig_transv_pos);
   fChain->SetBranchAddress("LorenzoFlag", &LorenzoFlag, &b_LorenzoFlag);
   fChain->SetBranchAddress("NpheCutSect", &NpheCutSect, &b_NpheCutSect);
   fChain->SetBranchAddress("trig_dcfid_flag", &trig_dcfid_flag, &b_trig_dcfid_flag);
   fChain->SetBranchAddress("MirrorCode", &MirrorCode, &b_MirrorCode);
   fChain->SetBranchAddress("nphot", &nphot, &b_nphot);
   fChain->SetBranchAddress("gam_Eec", gam_Eec, &b_gam_Eec);
   fChain->SetBranchAddress("gam_bet", gam_bet, &b_gam_bet);
   fChain->SetBranchAddress("gam_Pxtrue", gam_Pxtrue, &b_gam_Pxtrue);
   fChain->SetBranchAddress("gam_Pytrue", gam_Pytrue, &b_gam_Pytrue);
   fChain->SetBranchAddress("gam_Pztrue", gam_Pztrue, &b_gam_Pztrue);
   fChain->SetBranchAddress("gam_E", gam_E, &b_gam_E);
   fChain->SetBranchAddress("gam_Px", gam_Px, &b_gam_Px);
   fChain->SetBranchAddress("gam_Py", gam_Py, &b_gam_Py);
   fChain->SetBranchAddress("gam_Pz", gam_Pz, &b_gam_Pz);
   fChain->SetBranchAddress("gam_Etrue", gam_Etrue, &b_gam_Etrue);
   fChain->SetBranchAddress("gam_Etrue_Corr", gam_Etrue_Corr, &b_gam_Etrue_Corr);
   fChain->SetBranchAddress("gam_Pxtrue_Corr", gam_Pxtrue_Corr, &b_gam_Pxtrue_Corr);
   fChain->SetBranchAddress("gam_Pytrue_Corr", gam_Pytrue_Corr, &b_gam_Pytrue_Corr);
   fChain->SetBranchAddress("gam_Pztrue_Corr", gam_Pztrue_Corr, &b_gam_Pztrue_Corr);
   fChain->SetBranchAddress("gam_ecx", gam_ecx, &b_gam_ecx);
   fChain->SetBranchAddress("gam_ecy", gam_ecy, &b_gam_ecy);
   fChain->SetBranchAddress("gam_ecz", gam_ecz, &b_gam_ecz);
   fChain->SetBranchAddress("gam_ecu", gam_ecu, &b_gam_ecu);
   fChain->SetBranchAddress("gam_ecv", gam_ecv, &b_gam_ecv);
   fChain->SetBranchAddress("gam_ecw", gam_ecw, &b_gam_ecw);
   fChain->SetBranchAddress("gam_ect", gam_ect, &b_gam_ect);
   fChain->SetBranchAddress("gam_ecpath", gam_ecpath, &b_gam_ecpath);
   fChain->SetBranchAddress("gam_phi", gam_phi, &b_gam_phi);
   fChain->SetBranchAddress("gam_theta", gam_theta, &b_gam_theta);
   fChain->SetBranchAddress("npi0", &npi0, &b_npi0);
   fChain->SetBranchAddress("Pi0z_EVNT", Pi0z_EVNT, &b_Pi0z_EVNT);
   fChain->SetBranchAddress("Pi0mass_EVNT", Pi0mass_EVNT, &b_Pi0mass_EVNT);
   fChain->SetBranchAddress("opening_angle_evnt", opening_angle_evnt, &b_opening_angle_evnt);
   fChain->SetBranchAddress("Pi0energy_EVNT", Pi0energy_EVNT, &b_Pi0energy_EVNT);
   fChain->SetBranchAddress("Pi0z_ECPB", Pi0z_ECPB, &b_Pi0z_ECPB);
   fChain->SetBranchAddress("Pi0mass_ECPB", Pi0mass_ECPB, &b_Pi0mass_ECPB);
   fChain->SetBranchAddress("opening_angle_ecpb", opening_angle_ecpb, &b_opening_angle_ecpb);
   fChain->SetBranchAddress("Pi0energy_ECPB", Pi0energy_ECPB, &b_Pi0energy_ECPB);
   fChain->SetBranchAddress("Pi0mass_Corr", Pi0mass_Corr, &b_Pi0mass_Corr);
   fChain->SetBranchAddress("Pi0z_Corr", Pi0z_Corr, &b_Pi0z_Corr);
   fChain->SetBranchAddress("Pi0energy_Corr", Pi0energy_Corr, &b_Pi0energy_Corr);
   fChain->SetBranchAddress("Pi0pheta_Corr", Pi0pheta_Corr, &b_Pi0pheta_Corr);
   fChain->SetBranchAddress("Pi0phi_Corr", Pi0phi_Corr, &b_Pi0phi_Corr);
   fChain->SetBranchAddress("Pi0angle_gs", Pi0angle_gs, &b_Pi0angle_gs);
   fChain->SetBranchAddress("pT2_val_Corr", pT2_val_Corr, &b_pT2_val_Corr);
   fChain->SetBranchAddress("phi_val_Corr", phi_val_Corr, &b_phi_val_Corr);
   fChain->SetBranchAddress("pT1_val", pT1_val, &b_pT1_val);
   fChain->SetBranchAddress("opening_angle_Corr", opening_angle_Corr, &b_opening_angle_Corr);
   fChain->SetBranchAddress("G1", G1, &b_G1);
   fChain->SetBranchAddress("G2", G2, &b_G2);
   fChain->SetBranchAddress("eG1Angle", eG1Angle, &b_eG1Angle);
   fChain->SetBranchAddress("eG2Angle", eG2Angle, &b_eG2Angle);
   Notify();


}

Bool_t ec_gam::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef ec_gam_cxx

