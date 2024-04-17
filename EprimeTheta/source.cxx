#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TF1.h>
#include "ec_gam.h"

using namespace std;

int main(int argc, char *argv[])
{
TString target = "carbon";//set default to carbon`
int Ntoproc=0;
TChain *chain;
chain = new TChain("ec_gam","Analysis chain");
chain->Add("/lustre/expphy/volatile/clas/claseg2/taya/ID/ReducedID/carbon_id1.root");

int Ninchain = chain->GetEntries();
cout << "Number of events in chain : " << Ninchain;
if(Ntoproc==0)Ntoproc=Ninchain;
cout << " Will process : " << Ntoproc << endl;
ec_gam *selec = new ec_gam(0,Ntoproc,target);
chain->Process(selec,"",Ntoproc);
cout << "Back in main" << endl;

return 0;
}

