#define ec_gam_cxx
#include "ec_gam.h"
#include <TH2.h>
#include <TStyle.h>

void ec_gam::Begin(TTree *tree)
{
   TString option = GetOption();
   Init(tree);
   Nevta = 0;
   Nevt_found = 0;
   if(AvailableEvents==0)AvailableEvents = Int_t(tree->GetEntries());

   cout << "*******************************HELLO! ***********************************" << endl;
  
   cout << "Activating prefered styles... "; cout.flush();
   gROOT->SetStyle("Plain");
   gStyle->SetPalette(1);
  // gStyle->SetOptFit(1111);
  // gStyle->SetOptStat(11);
  // gStyle->SetOptTitle(1);
   cout << "Done !" << endl << "Allocations/Initializations"; cout.flush();
   cout << endl << "The time is now ";cout.flush(); gSystem->Exec("date");cout<<endl;
  
   NewFile = new TFile("NewFile.root","recreate");
   TreeGam = new TTree("ec_gam","Tree with electron variables");
   MakeBranch();
   cout<< "Branch made" <<endl;  

   BookHists();

   cout << "Done !" << endl;
   cout << "looping" << endl;
}

void ec_gam::SlaveBegin(TTree *tree)
{
   TString option = GetOption();
   Init(tree);
}

Bool_t ec_gam::Process(Long64_t entry)
{
   fChain->GetTree()->GetEntry(entry);
   TString option = GetOption();
   Nevta++;
  
   Int_t evtsPrint = Int_t(TMath::Floor(Nevta/(AvailableEvents/100)));
   evtsPrint = evtsPrint*Int_t(AvailableEvents/100);
   if(Nevta==evtsPrint){
   cout << "\r"<<TMath::FloorNint(Nevta/(AvailableEvents/100.))<<"\%";cout.flush();
   }
   
   Calculate();
   FillTables();
   TreeGam->Fill();

   return kTRUE;
}

void ec_gam::SlaveTerminate()
{

}

void ec_gam::Terminate()
{

  WriteHists();
  cout<<" hist written " << endl;
  WriteTxt();
  cout<< " Txt is written " << endl;
  NewFile->Close();
  cout<<"file is closed"<<endl;
  

  cout << "Over. Ciao, baby!" << endl;
}


void ec_gam::Calculate(void){
// here calculate Eprime and Theta
 Eprime = Ebeam - trig_nu;
 Theta = 2*(TMath::ASin(TMath::Sqrt(trig_Q2/(4*Ebeam*Eprime)))) ;  
 Theta = TMath::RadToDeg()*Theta;

 //cout<< "trig_nu = "<< trig_nu<< " Eprime=" << Eprime << " Theta=" << Theta<<endl;

}

void ec_gam::BookHists(void){
 qbin = -1;
 vbin = -1;
 h_etheta = new TH2F("h_etheta","Theta_Eprime",200,0.5,3.,200,15,55.); 
 for(int  nq=0;  nq<nqbins; nq++){
  for(int  nv=0; nv<nvbins;nv++){
   h_etheta_b[nq][nv] = new TH2F(Form("h_etheta_b%d%d",nq,nv),Form("h_etheta_b%d%d",nq,nv),200,0.5,3.,200,15,55.);
   }
  }

}


void ec_gam::FillTables(void){
 qbin = Get_q_bin(trig_Q2);
 vbin =  Get_v_bin(trig_nu); 
 if(qbin>-1&&qbin<nqbins &&vbin>-1&&vbin<nvbins){
  h_etheta->Fill(Eprime, Theta);
  h_etheta_b[qbin][vbin]->Fill(Eprime, Theta);
 }

}

void ec_gam::WriteHists(void){

 gStyle->SetOptTitle(0);
 TCanvas *c = new TCanvas("c","c",700,500);
 h_etheta->Draw("colz");
 h_etheta->GetXaxis()->SetTitle("E'(GeV)");
 h_etheta->GetYaxis()->SetTitle("#theta (deg)");


//////////////// nu centers //////////////
/*
 for(int nv=0; nv<nvbins;nv++){
 TLine *nu = new TLine(Ebeam -center_v[nv],15., Ebeam - center_v[nv],55);
 if(nv<2) nu->SetLineColor(nv+1);
 else nu->SetLineColor(kGray+1);
 nu->SetLineWidth(2); 
 nu->Draw("same");
 }
*/
//////////////// nu borders //////////////

 for(int nv=0; nv<nvbins+1;nv++){
 TLine *nu = new TLine(Ebeam - v_values[nv],15., Ebeam - v_values[nv],55);
 if(nv==0) nu->SetLineColor(1);
 else nu->SetLineColor(kGray-nv+4);
 nu->SetLineWidth(2); 
 nu->Draw("same");
 }

///////////// Q2 lines /////////////////
/*
 for(int nq=0; nq<nqbins;nq++){
  TF1 *Q2lines = new TF1(Form("Q2_%d",nq+1),"TMath::RadToDeg()*2*TMath::ASin(TMath::Sqrt([0]/(4*[1]*x)))", 0.7, 2.8);
  Q2lines->SetParameters(center_q[nq], Ebeam);
  if(nq<2)Q2lines->SetLineColor(nq+1);
  else Q2lines->SetLineColor(kGray+1);
  Q2lines->SetLineWidth(2);
  Q2lines->Draw("same");  
 }
*/
///////////// Q2 lines /////////////////
 for(int nq=0; nq<nqbins+1;nq++){
  TF1 *Q2lines = new TF1(Form("Q2_%d",nq+1),"TMath::RadToDeg()*2*TMath::ASin(TMath::Sqrt([0]/(4*[1]*x)))", 0.764, 2.814);
  Q2lines->SetParameters(q_values[nq], Ebeam);
  if(nq==0)Q2lines->SetLineColor(1);
  Q2lines->SetLineColor(kGray-nq+4);
  Q2lines->SetLineWidth(2);
  Q2lines->Draw("same");  
 }

///////////////////////////////////////////

 c->Print("~/secure/ExternalsRC/ThetaEprime_Q2NuBorders.png");
 c->Clear();

///////////////////////// Binned //////////////////////////
 gStyle->SetOptTitle(11111);
 TCanvas *c1 = new TCanvas("c1","c1",400*3,300*3);
 c1->Divide(3,3, 0.0001,0.0001);
 for(int  nq=0;  nq<nqbins; nq++){
   for(int  nv=0; nv<nvbins;nv++){
    c1->cd(nq+1+3*nv);
    h_etheta_b[nq][nv]->Draw("colz");
    h_etheta_b[nq][nv]->GetXaxis()->SetTitle("E'(GeV)");
    h_etheta_b[nq][nv]->GetYaxis()->SetTitle("#theta(deg)");
    
    mean_theta[nq][nv] = h_etheta_b[nq][nv]->GetMean(2);
    mean_Eprime[nq][nv] =h_etheta_b[nq][nv]-> GetMean(1);    
     
    cout << nq << " " << nv << " <Q2>="<< center_q[nq] << " <v>="<<center_v[nv] << " <theta>="<< h_etheta_b[nq][nv]->GetMean(2)
                            << " <E'>="<< h_etheta_b[nq][nv]-> GetMean(1)<<endl;
    }
   } 
 c1->Print("~/secure/ExternalsRC/ThetaEprime_Q2NuBorders_binned.png");
 c1->Clear();
}
////////////////// TREE /////////////////////////////////

void ec_gam::MakeBranch(void){
 if(TreeGam!=NULL){
  if(qbin>-1&&qbin<nqbins &&vbin>-1&&vbin<nvbins){
  TreeGam->Branch("trig_theta",&trig_theta,"trig_theta/F");
  TreeGam->Branch("trig_Q2",&trig_Q2,"trig_Q2/F");
  TreeGam->Branch("trig_nu",&trig_nu,"trig_nu/F");
  TreeGam->Branch("trig_mom",&trig_mom,"trig_mom/F");
  TreeGam->Branch("Theta", &Theta,"Theta/F");
  TreeGam->Branch("Eprime", &Eprime,"Eprime/F");
  }
 }
 
 return;
}

void ec_gam::WriteTxt(void){
ofstream outfile;
outfile.open("qv_Eprime_Theta.txt");
for(int  nq=0;  nq<nqbins; nq++){
  for(int  nv=0; nv<nvbins;nv++){
     //outfile << nq << " " << nv<< " " << mean_theta[nq][nv]<< " "<<mean_Eprime[nq][nv]<<endl;
     outfile << mean_theta[nq][nv]<< " "<<mean_Eprime[nq][nv]<<endl;
   }
  }
outfile.close();

} // writeTXT

/////////////////////////////////////////////////////////
int  ec_gam::Get_q_bin(float q){
  if(q<qmin||q>qmax)return -1;
   for(int nq=0; nq<nqbins; nq++){
    if(q>q_values[nq]&&q<q_values[nq+1]) return nq;
		      }
   return -1;
}

 int  ec_gam::Get_v_bin(float v){
  if(v<vmin||v>vmax)return -1;
   for(int nv=0; nv<nvbins; nv++){
    if(v>v_values[nv]&&v<v_values[nv+1]) return nv;
		       }
   return -1;
 }
/*
 Double_t ec_gam::thetafnc(Double_t *x, Double_t *par)
{
   Float_t xx =x[0];
   Double_t f = TMath::ASin(TMath::Sqrt(par[0]/(5.014*xx))); 
   return f;
   //return TMath::RadToDeg()*f;
}
*/




