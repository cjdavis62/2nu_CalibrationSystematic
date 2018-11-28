
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "TROOT.h"
#include "TTree.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TCut.h"
#include <iostream>
#include "THStack.h"
#include "TMath.h"

using namespace std;

// Read in the data file with all the shifts at certain energies
// Take the shift as being symmetric (conservative)
double ReadShift() {
  stringstream ss;
  string line;
  
  double energy, energyShift, energyShiftSigma;
  double max_shift;
  
  ifstream shiftFile("EnergyShift.txt");
  while (getline(shiftFile, line)) {
    ss.str(line);
    ss >> energy >> energyShift >> energyShiftSigma;
    ss.str("");
    ss.clear();
    cout << energy << " " << energyShift << " " << energyShiftSigma << endl;

    max_shift = TMath::Max(TMath::Abs(energyShift + energyShiftSigma), TMath::Abs(energyShift - energyShiftSigma));
    cout << max_shift << endl;
  }
    return max_shift;
}


// Given a certain energy to the events, give a keV shift
double EnergyEdit(double energy_data) {

  
}

void AddCalibrationFunction() {

  ReadShift();
  
  // Read list of filenames
  
  // Find the MC file to read
  TFile * mcFile = new TFile("/projects/cuore/simulation/MC_test.root", "update");
  TTree * mcTree = (TTree*) mcFile->Get("outTree");

  // Tell the reader which branches are interesting
  mcTree->SetBranchStatus("*", 0);
  mcTree->SetBranchStatus("Ener2", 1);
  mcTree->SetBranchStatus("ESum2", 1);
  mcTree->SetBranchStatus("Multiplicity", 1);
  mcTree->SetBranchStatus("MultipletIndex", 1);
  
  // Get the branch addresses
  double Ener2_data;
  double ESum2_data;
  Short_t Multiplicity;
  Short_t MultipletIndex;

  mcTree->SetBranchAddress("Ener2", &Ener2_data);
  mcTree->SetBranchAddress("ESum2", &ESum2_data);
  mcTree->SetBranchAddress("Multiplicity", &Multiplicity);
  mcTree->SetBranchAddress("MultipletIndex", &MultipletIndex);

  // Create new branches for the MC file

  
  // Grab all the events in the file
  Long64_t e = 0;
  while(e < mcTree->GetEntries()) {
    mcTree->GetEntry(e);

    // Treat high multiplicity events seperately
    if(Multiplicity == 1) {
      // Edit the energies for totalEnergy and Energy
      // Fill the tree
    }
    // If Multiplicity > 1, need to look at all the events now to adjust ESum
    else {
      for (int i = 1; i <= Multiplicity; i++) {
	// Edit the energies for totalEnergy and Energy of each event
      }
      // Sum the new totalEnergy for the events
      // Fill the tree for each
      
    }
    //cout << e << endl;
  }
  
  

  // Rewrite the MC file


}
