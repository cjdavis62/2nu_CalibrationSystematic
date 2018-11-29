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
#include <sstream>
#include <fstream>


using namespace std;

// Read in the data file with all the shifts at certain energies
// Take the shift as being symmetric (conservative)
std::pair<std::vector<double>, std::vector<double> > ReadShift(std::string Shifted_file_name) {
  // structures to read in the data
  stringstream ss;
  string line;

  double energy, energyShift, energyShiftSigma;
  double max_shift;

  // store all the data in vectors
  std::vector<double> energy_vector, energyShift_vector;

  // Read the file
  ifstream shiftFile( const_cast<char*>(Shifted_file_name.c_str()));
  
  while (getline(shiftFile, line)) {
    ss.str(line);
    ss >> energy >> energyShift >> energyShiftSigma;
    ss.str("");
    ss.clear();

    //cout << energy << " " << energyShift << " " << energyShiftSigma << endl;

    // use the biggest delta as the max_shift here. Make it symmetric to be conservative
    max_shift = TMath::Max(TMath::Abs(energyShift + energyShiftSigma), TMath::Abs(energyShift - energyShiftSigma));
    //cout << max_shift << endl;
    energy_vector.push_back(energy);
    energyShift_vector.push_back(max_shift);
  }
  // return the two vectors as a pair
  return std::make_pair(energy_vector, energyShift_vector);
}

std::vector<std::string> ReadFilenames(std::string Filename){
  stringstream ss;
  std::string line;
  std::string MC_file;
  std::vector<std::string> file_vector;
  
  ifstream file(const_cast<char*>(Filename.c_str()));
  while(getline(file, line)) {
    ss.str(line);
    ss >> MC_file;
    
    file_vector.push_back(MC_file);
  }
  return file_vector;
}

// Given a certain energy to the events, give a keV shift
double EnergyEdit(double energy_data, int up_or_down) {

  return 0;
  
}

void AddCalibrationFunction() {

  // The file containing the energy shifts
  std::string shifted_filename = "EnergyShift.txt";

  // Read the file and record the data as a vector in a pair
  // Format looks like
  // || energy | shifted energy ||
  std::pair<std::vector<double>, std::vector<double> > shifted_pair = ReadShift(shifted_filename);

  std::vector<double> energy_vector = shifted_pair.first;
  std::vector<double> energyShift_vector = shifted_pair.second;

  cout << energy_vector[0] << " " << energyShift_vector[0] << endl;
  
  // Read list of filenames
  // std::vector<std::string> ReadFilenames(std::string Filename)
  
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
  double Ener2_data, ESum2_data;
  Short_t Multiplicity;
  Short_t MultipletIndex;
  
  mcTree->SetBranchAddress("Ener2", &Ener2_data);
  mcTree->SetBranchAddress("ESum2", &ESum2_data);
  mcTree->SetBranchAddress("Multiplicity", &Multiplicity);
  mcTree->SetBranchAddress("MultipletIndex", &MultipletIndex);

  // Create new tree for the MC file
  double Ener2_shift_up, Ener2_shift_down, ESum2_shift_up, ESum2_shift_down;
  TTree* friendTree = new TTree("shiftedEnergyTree", "");
  friendTree->Branch("Ener2_shift_up", &Ener2_shift_up, "Ener2_shift_up/D");
  friendTree->Branch("Ener2_shift_down", &Ener2_shift_down, "Ener2_shift_down/D");
  friendTree->Branch("ESum2_shift_up", &ESum2_shift_up, "ESum2_shift_up/D");
  friendTree->Branch("ESum2_shift_down", &ESum2_shift_down, "ESum2_shift_down/D");
  
  // The new energies
  double ESum_shift_up, ESum_shift_down;
  vector<double> energy_shift_up, energy_shift_down;
  
  bool previously_edited = false;
  
  // Grab all the events in the file
  Long64_t e = 0;
  while (e < mcTree->GetEntries()) {

    // Start the loop for each new set of Multiplets
    if (energy_shift_up.empty())
      {
	// Set the Total Energy to 0
	ESum_shift_up = 0;
	ESum_shift_down = 0;
	int i = 0;
	while (true) {
	  // Get the entry in the root file
	  mcTree->GetEntry(e + i);
	  // Edit the energies into a vector. The vector will have size == Multiplicity by the end of each set
	  energy_shift_up.push_back(EnergyEdit(Ener2_data, 1));
	  energy_shift_down.push_back(EnergyEdit(Ener2_data, -1));
	  // Add to the sum
	  ESum_shift_up += energy_shift_up.back();
	  ESum_shift_down -= energy_shift_down.back();
	  i++;
	  if (i == Multiplicity) break;
	}
      }

    // Pop out the first element of the vectors to get the energy of the event
    Ener2_shift_up = energy_shift_up[0];
    Ener2_shift_down = energy_shift_down[0];
    energy_shift_up.erase(energy_shift_up.begin());
    energy_shift_down.erase(energy_shift_down.begin());

    // Fill the tree for each
    friendTree->Fill();
    // Go to next event
    e++;
  }

  // Rewrite the MC file
  friendTree->Write();
  mcFile->Close();
  delete mcFile;
}
