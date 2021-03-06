/* this script reads in a list of MC files and applies a positive and negative calibration shift to them
compile with:
> g++ -o AddCalibrationFunction AddCalibrationFunction.cc `root-config --cflags --glibs`
and run with
> ./AddCalibrationFunction mc_file

The text file containing the data for the energy shift should be in the same directory and named 'EnergyShift.txt'
*/

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
#include <iomanip>

using namespace std;

// Read in the data file with all the shifts at certain energies
// Returns a pair of maps, one with the energy shifted up, the other down
std::pair<std::map<double, double>, std::map<double, double> > ReadShift(std::string Shifted_file_name) {
  // structures to read in the data
  stringstream ss;
  string line;

  double energy, energyShift, energyShiftSigma;
  double shift_up;
  double shift_down;

  // Do instead as a map
  std::map<double, double> energy_shift_up;
  std::map<double, double> energy_shift_down;

  // Read the file
  ifstream shiftFile( const_cast<char*>(Shifted_file_name.c_str()));
  
  while (getline(shiftFile, line)) {
    ss.str(line);
    ss >> energy >> energyShift >> energyShiftSigma;
    ss.str("");
    ss.clear();
        
    // use the biggest delta as the max_shift here. Make it symmetric to be conservative
    shift_up = energyShift + energyShiftSigma;
    shift_down = energyShift - energyShiftSigma;
    
    energy_shift_up[energy] = shift_up;
    energy_shift_down[energy] = shift_down;
  }
  shiftFile.close();
  // return the map of energies as a pair
  return make_pair(energy_shift_up, energy_shift_down);
}

// Read the txt file storing the name and path to the MC files. In the text file, the first line is the path, and the lines after are each of the files
// Returns the full path to the MC files as a vector
std::vector<std::string> ReadFilenames(std::string Filename) {
  stringstream ss;
  std::string line;
  std::string MC_file;
  std::string MC_dir;
  std::vector<std::string> file_vector;
  
  ifstream file(const_cast<char*>(Filename.c_str()));
  // Read in the first line to get the directory
  int i = 0;
  while(getline(file, line)) {
    if (i == 0)
      {
	ss.str(line);
	ss >> MC_dir;
	i++;
      }
    else
      {
	ss.str(line);
	ss >> MC_file;
	file_vector.push_back(MC_dir + "/" + MC_file);
      }
    ss.str("");
    ss.clear();
  }
  file.close();
  return file_vector;
}

// Given a certain energy to the events, give a keV shift
double EnergyEdit(std::map<double, double> energy_shift_map, double energy_data) {

  // Because the energy shift map has the data as 0.5, 1.5, etc, and floor and ceil functions are great, shift all the data by 0.5 up so we can use them
  energy_data = energy_data + 0.5;

  double energy_data_low = std::floor(energy_data) - 0.5;
  double energy_data_high = std::ceil(energy_data) - 0.5;

  // Restore the data to its original value
  energy_data = energy_data - 0.5;

  // linearly interpolate the values of the function between these two
  double energy_shift_interpolated = energy_shift_map[energy_data_low] + (energy_shift_map[energy_data_high] - energy_shift_map[energy_data_low]) * (energy_data - energy_data_low);
  
  return energy_data + energy_shift_interpolated;
}
// Adds the shift to the data file
void AddCalibrationFunction(std::string MC_file, std::map<double, double> energy_shift_up_map, std::map<double, double> energy_shift_down_map) {

  std::cout << "Reading file: " << MC_file << std::endl;
  
  // Get the MC file to read
  TFile * mcFile = new TFile(const_cast<char*>(MC_file.c_str()), "update");
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
  double Ener2_shift_up, Ener2_shift_down;
  double ESum2_shift_up = 0;
  double ESum2_shift_down = 0 ;

  // Make sure the TTree doesn't have any friends with the same name. It can only be friends with the new kid.
  TTree* friendTree;
  if ((friendTree = (TTree*) mcFile->Get("shiftedEnergyTree")) != NULL) {
    mcTree->RemoveFriend(friendTree);
    friendTree->Delete("all");
  }
  // Make the new friend TTree
  friendTree = new TTree("shiftedEnergyTree", "");
  friendTree->Branch("Ener2_shift_up", &Ener2_shift_up, "Ener2_shift_up/D");
  friendTree->Branch("Ener2_shift_down", &Ener2_shift_down, "Ener2_shift_down/D");
  friendTree->Branch("ESum2_shift_up", &ESum2_shift_up, "ESum2_shift_up/D");
  friendTree->Branch("ESum2_shift_down", &ESum2_shift_down, "ESum2_shift_down/D");
  
  // Store the energies in a vector so that we can deal with Multiplicity > 1 events
  vector<double> energy_shift_up, energy_shift_down;
  
  bool previously_edited = false;
  
  // Grab all the events in the file
  Long64_t e = 0;
  double percent;
  int i;
  while (e < mcTree->GetEntries()) {

    // Tell the user how far we are through the MC file every 10k events
    percent = (double(e) / mcTree->GetEntries()) * 100.0;
    if (e % 10000 == 0) {
      std::cout << "Entries read: " << e << " of " << mcTree->GetEntries() << " | " << std::fixed <<  std::setprecision(4) << percent << "% complete" << "\r" << std::flush;
    }
    
    // Start the loop for each new set of Multiplets
    if (energy_shift_up.empty())
      {
	// Set the Total Energy to 0
	ESum2_shift_up = 0;
	ESum2_shift_down = 0;
	i = 0;
	while (true) {
	  //std::cout << "i: " << i << " Multiplicity: " << int(Multiplicity) << std::endl;
	  // Get the entry in the root file
	  mcTree->GetEntry(e + i);
	  // Edit the energies into a vector. The vector will have size == Multiplicity by the end of each set
	  energy_shift_up.push_back(EnergyEdit(energy_shift_up_map, Ener2_data));
	  energy_shift_down.push_back(EnergyEdit(energy_shift_down_map, Ener2_data));
	  // Add to the sum
	  ESum2_shift_up += energy_shift_up.back();
	  ESum2_shift_down += energy_shift_down.back();
	  // Check if we need to stop or continue
	  if (i + 1 >= int(Multiplicity))
	    {
	      break;
	    }
	  else
	    {
	      i++;
	    }
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
  std::cout << "\n";
  
  // Rewrite the MC file
  friendTree->Write();
  mcFile->Close();
  delete mcFile;
}


int main(int argc, char **argv) {
  
  if (argc !=2) {
    std::cout << "usage: ./AddCalibrationFunction MC_file.root" << std::endl;
    exit(1);
  }
  
  std::string mc_file = argv[1]; // The mc file is passed through the arguments of the function by the python script
  
  // The file containing the energy shifts
  std::string shifted_filename = "EnergyShift.txt";

  // Read the file and record the data as a map
  // Format looks like
  // || energy | shifted energy ||
  std::pair<std::map<double, double>, std::map<double, double> >  energy_shift_map_pair = ReadShift(shifted_filename);
  std::map<double, double> energy_shift_up = energy_shift_map_pair.first;
  std::map<double, double> energy_shift_down = energy_shift_map_pair.second;
  
  AddCalibrationFunction(mc_file, energy_shift_up, energy_shift_down);
  
}
