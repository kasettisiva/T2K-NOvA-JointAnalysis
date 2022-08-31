#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TVectorD.h"

#include <iostream>
#include <vector>
#include <string>
using namespace std;

const int ToyXpOscThrowOffset = 1; // we need this offset because for asimov fits ToyXpOscThrow will be -1

void MergeToyMargTemplates(const char *inpattern, const char *outname) {
TChain *intree = new TChain("MargTemplate");
intree->Add(inpattern);

Double_t delta;
Double_t mh;
Int_t Index;
Int_t ToyXpOscThrow;
Double_t ToyXpWeight = 1.;
Double_t AvNLLtot;
Double_t AvNLLtot2;
Double_t AvNLLtot_NOvA;
Double_t AvNLLtot_T2KNOvA;

vector<string> oscParamNames;
oscParamNames.push_back("sin2213");
oscParamNames.push_back("delta");
oscParamNames.push_back("dm2");
oscParamNames.push_back("sin223");
oscParamNames.push_back("dm2_bar");
oscParamNames.push_back("sin223_bar");
oscParamNames.push_back("mh");
oscParamNames.push_back("smear_fact");

int nosc = oscParamNames.size();
vector<bool> useOscParam(nosc);
vector<double> oscParams(nosc);
vector<vector<double> > oscParamAtIndex(nosc); // [param][index]
int first_used_io = -1;
for (int io = 0; io < nosc; io++) {
    if (intree->GetBranch(oscParamNames.at(io).c_str()) != NULL) {
        cout << "Tree has " << oscParamNames.at(io) << endl;
        useOscParam.at(io) = true;
        intree->SetBranchAddress(oscParamNames.at(io).c_str(), &oscParams.at(io));

        if (first_used_io < 0) {
          first_used_io = io;
        }
    }
}

intree->SetBranchAddress("Index", &Index);
intree->SetBranchAddress("ToyXpOscThrow", &ToyXpOscThrow);
if (intree->GetBranch("ToyXpWeight") != NULL) {
    intree->SetBranchAddress("ToyXpWeight", &ToyXpWeight);
}
intree->SetBranchAddress("AvNLLtot", &AvNLLtot);
bool has_AvNLLtot2 = (intree->GetBranch("AvNLLtot2") != NULL);
if (has_AvNLLtot2) {
    intree->SetBranchAddress("AvNLLtot2", &AvNLLtot2);
}
intree->SetBranchAddress("AvNLLtot_NOvA", &AvNLLtot_NOvA);
intree->SetBranchAddress("AvNLLtot_T2KNOvA", &AvNLLtot_T2KNOvA);

Long_t Nentries = intree->GetEntries();
if (Nentries == 0) {
    cerr << "Zero entries." << endl;
}

std::cout << " Nentries " << Nentries << std::endl;

int Nbins = -1;
for (int ev = 0; ev < Nentries; ev++) {
  intree->GetEntry(ev);
  
  if (Index != (int)oscParamAtIndex.at(first_used_io).size()) {
    // Index is always incremental, so we can find the max
    // by detecting when it becomes smaller
    
    // also this check makes sure that the deltaAtIndex key
    // matches the Index value
    Nbins = oscParamAtIndex.at(first_used_io).size();
    printf("Number of gridpoints is %u\n", Nbins);
    break;
  }
  
  for (int io = 0; io < nosc; io++) {
    if (useOscParam.at(io)) {
      oscParamAtIndex.at(io).push_back(oscParams.at(io));
    }
  }
}

if (first_used_io < 0) {
  cout << "No grid?" << endl;
  Nbins = 1;
}

std::vector<double> emptyGrid(Nbins);
std::vector<double> likelihoodOffset(Nbins);
std::vector<bool>   didSetLikelihoodOffset(Nbins);
std::vector<double> likelihoodOffset_NOvA(Nbins);
std::vector<double> likelihoodOffset_T2KNOvA(Nbins);

// this ToyXpIndex is something we just introduce here for the vector indexing,
// and is unrelated to the ToyXp value used in the MargTemplate tree
std::vector<int> ToyXpIndexForToyXpOscThrow;
std::vector<int> ToyXpOscThrowForToyXpIndex;

std::vector<std::vector<double> > toyEntries;    // [toy][index]
std::vector<std::vector<double> > toyLikelihood; // [toy][index]
std::vector<std::vector<double> > toyLikelihood2; // [toy][index]
std::vector<std::vector<double> > toyLikelihood_NOvA; // [toy][index]
std::vector<std::vector<double> > toyLikelihood_T2KNOvA; // [toy][index]
std::vector<double> toyWeight; // [toy]

const int kUndefinedToyXpIndex = -1;

for (int ev = 0; ev < Nentries; ev++) {
  intree->GetEntry(ev);
  //std::cout << " Index " << Index << std::endl;
  if (ev % 1000000 == 0) { printf("%ld/%ld\n", ev, Nentries); }
  
  if (ToyXpOscThrowOffset+ToyXpOscThrow >= (int)ToyXpIndexForToyXpOscThrow.size()) {
    ToyXpIndexForToyXpOscThrow.resize(ToyXpOscThrowOffset+ToyXpOscThrow+1, kUndefinedToyXpIndex);
  }
  
  if (ToyXpIndexForToyXpOscThrow.at(ToyXpOscThrowOffset+ToyXpOscThrow) == kUndefinedToyXpIndex) {
    // need to add new entry
    int ToyXpIndex = ToyXpOscThrowForToyXpIndex.size();
    ToyXpIndexForToyXpOscThrow.at(ToyXpOscThrowOffset+ToyXpOscThrow) = ToyXpIndex;
    ToyXpOscThrowForToyXpIndex.push_back(ToyXpOscThrow);
    // also enlarge the matrices
    toyEntries   .push_back(emptyGrid);
    toyLikelihood.push_back(emptyGrid);
    toyLikelihood_NOvA.push_back(emptyGrid);
    toyLikelihood_T2KNOvA.push_back(emptyGrid);
    toyLikelihood2.push_back(emptyGrid);
    toyWeight    .push_back(ToyXpWeight);
  }
  
  if (!didSetLikelihoodOffset.at(Index)) {
    //std::cout << " AvNLLtot " << AvNLLtot << std::endl;
    likelihoodOffset.at(Index) = AvNLLtot;
    likelihoodOffset_NOvA.at(Index) = AvNLLtot_NOvA;
    likelihoodOffset_T2KNOvA.at(Index) = AvNLLtot_T2KNOvA;
    didSetLikelihoodOffset.at(Index) = true;
  }
  
  int ToyXpIndex = ToyXpIndexForToyXpOscThrow.at(ToyXpOscThrowOffset+ToyXpOscThrow);
  toyEntries   .at(ToyXpIndex).at(Index)++;
  toyLikelihood.at(ToyXpIndex).at(Index) += exp(-(AvNLLtot-likelihoodOffset.at(Index)));
  toyLikelihood_NOvA.at(ToyXpIndex).at(Index) += exp(-(AvNLLtot_NOvA-likelihoodOffset_NOvA.at(Index)));
  toyLikelihood_T2KNOvA.at(ToyXpIndex).at(Index) += exp(-(AvNLLtot_T2KNOvA-likelihoodOffset_T2KNOvA.at(Index)));
  if (has_AvNLLtot2) {
      toyLikelihood2.at(ToyXpIndex).at(Index) += exp(-(AvNLLtot2-2.*likelihoodOffset.at(Index)));
  }
}

std::cout << " toyEntries " << toyEntries.size() << "  " << toyEntries[0].size()<< std::endl;

std::vector<double> toyAvNLLtot(Nbins);
std::vector<double> toyAvNLLtot2(Nbins);
std::vector<double> toyAvNLLtot_NOvA(Nbins);
std::vector<double> toyAvNLLtot_T2KNOvA(Nbins);

TFile *outfile = new TFile(outname, "RECREATE");
TTree *outtree = new TTree("MargTemplatePerToy", "");
int MinNLLIndex; // index of NLL minimum
outtree->Branch("ToyXpOscThrow", &ToyXpOscThrow, "ToyXpOscThrow/I");
outtree->Branch("MinNLLIndex", &MinNLLIndex, "MinNLLIndex/I");
outtree->Branch("ToyXpWeight", &ToyXpWeight, "ToyXpWeight/D");
outtree->Branch("AvNLLtot", &toyAvNLLtot);
outtree->Branch("AvNLLtot_NOvA", &toyAvNLLtot_NOvA);
outtree->Branch("AvNLLtot_T2KNOvA", &toyAvNLLtot_T2KNOvA);
if (has_AvNLLtot2) {
    outtree->Branch("AvNLLtot2", &toyAvNLLtot2);
}

TTree *legacyTree = new TTree("MargTemplate", ""); // this is the same information in MargTemplate format for plotting tools
Int_t ToyXp = -1;
legacyTree->Branch("Index", &Index, "Index/I");
legacyTree->Branch("ToyXp", &ToyXp, "ToyXp/I");
legacyTree->Branch("ToyXpOscThrow", &ToyXpOscThrow, "ToyXpOscThrow/I");
legacyTree->Branch("ToyXpWeight", &ToyXpWeight, "ToyXpWeight/D");
legacyTree->Branch("AvNLLtot", &AvNLLtot, "AvNLLtot/D");
if (has_AvNLLtot2) {
    legacyTree->Branch("AvNLLtot2", &AvNLLtot2, "AvNLLtot2/D");
}
legacyTree->Branch("AvNLLtot_NOvA", &AvNLLtot_NOvA, "AvNLLtot_NOvA/D");
legacyTree->Branch("AvNLLtot_T2KNOvA", &AvNLLtot_T2KNOvA, "AvNLLtot_T2KNOvA/D");
for (int io = 0; io < nosc; io++) {
    if (useOscParam.at(io)) {
        legacyTree->Branch(oscParamNames.at(io).c_str(), &oscParams.at(io), Form("%s/D", oscParamNames.at(io).c_str()));
    }
}

int previousN = -1;
for (int ToyXpIndex = 0; ToyXpIndex < (int) toyEntries.size(); ToyXpIndex++) {
  ToyXp = ToyXpIndex; // for legacyTree
  ToyXpOscThrow = ToyXpOscThrowForToyXpIndex.at(ToyXpIndex);
  if (ToyXpIndex < 2) {
    printf("ToyXpOscThrow = %d\n", ToyXpOscThrow);
  }
  
  for (int ib = 0; ib < Nbins; ib++) {
    Index = ib;
    AvNLLtot = 0.;
    AvNLLtot_NOvA = 0.;
    AvNLLtot_T2KNOvA = 0.;
    AvNLLtot2 = 0.;
    double N = toyEntries.at(ToyXpIndex).at(ib);
    if (previousN != -1) {
      if (previousN != N) { cerr << "Inconsistent grid sizes?" << endl; return; }
    }
    previousN = N;
    if (toyEntries.at(ToyXpIndex).at(ib) > 0) {
      double L = toyLikelihood.at(ToyXpIndex).at(ib);
      double L_NOvA = toyLikelihood_NOvA.at(ToyXpIndex).at(ib);
      double L_T2KNOvA = toyLikelihood_T2KNOvA.at(ToyXpIndex).at(ib);

      AvNLLtot  = likelihoodOffset.at(ib) - log(L/N);
      AvNLLtot_NOvA  = likelihoodOffset_NOvA.at(ib) - log(L_NOvA/N);
      AvNLLtot_T2KNOvA  = likelihoodOffset_T2KNOvA.at(ib) - log(L_T2KNOvA/N);

      toyLikelihood.at(ToyXpIndex).at(ib) = AvNLLtot;
      toyLikelihood_NOvA.at(ToyXpIndex).at(ib) = AvNLLtot_NOvA;
      toyLikelihood_T2KNOvA.at(ToyXpIndex).at(ib) = AvNLLtot_T2KNOvA;

      if (has_AvNLLtot2) {
          double L2 = toyLikelihood2.at(ToyXpIndex).at(ib);
          AvNLLtot2 = 2.*likelihoodOffset.at(ib) - log(L2/N);
          toyLikelihood2.at(ToyXpIndex).at(ib) = AvNLLtot2;
      }
    }
    for (int io = 0; io < nosc; io++) {
      if (useOscParam.at(io)) {
        oscParams.at(io) = oscParamAtIndex.at(io).at(ib);
      }
    }
    legacyTree->Fill();

    if (ToyXpIndex < 2 && ib < 3) {
      printf("  llh[%d] = %g\n", ib, toyLikelihood.at(ToyXpIndex).at(ib));
    }
  }
  
  ToyXpWeight = toyWeight.at(ToyXpIndex);
  toyAvNLLtot = toyLikelihood.at(ToyXpIndex);
  toyAvNLLtot_NOvA = toyLikelihood_NOvA.at(ToyXpIndex);
  toyAvNLLtot_T2KNOvA = toyLikelihood_T2KNOvA.at(ToyXpIndex);

  toyAvNLLtot2 = toyLikelihood2.at(ToyXpIndex);
  MinNLLIndex = TMath::LocMin(Nbins, toyAvNLLtot.data());
  outtree->Fill();
}

// let's also make one MargTemplate 


for (int io = 0; io < nosc; io++) {
  if (useOscParam.at(io)) {
    TVectorD v_oscParamAtIndex(Nbins, oscParamAtIndex.at(io).data());
    v_oscParamAtIndex.Write(Form("%sGrid", oscParamNames.at(io).c_str()));
  }
}
outtree->Write();
legacyTree->Write();
outfile->Close();
cout << "Wrote to " << outfile->GetName() << endl;
}
