/* this script reads in a list of MC files and applies a positive and negative calibration shift to them
compile with:
> g++ -o AddCalibrationFunction_newfile AddCalibrationFunction_newfile.cc `root-config --cflags --glibs`
and run with
> ./AddCalibrationFunction_newfile mc_file split number_of_splits

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
#include <string>

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
void AddCalibrationFunction(std::string MC_file, std::string directory, std::map<double, double> energy_shift_up_map, std::map<double, double> energy_shift_down_map, int split, int number_of_splits) {

  std::cout << "Reading file: " << MC_file << std::endl;
  
  // Get the MC file to read
  TFile * mcFile = new TFile(const_cast<char*>(MC_file.c_str()), "READ");
  std::string outfile_name;
  stringstream ss;
  ss << split;
  outfile_name = directory + "/" + "root_" +  ss.str();
  TFile * outFile = new TFile(const_cast<char*>(outfile_name.c_str()), "RECREATE");
  
  TTree * mcTree = (TTree*) mcFile->Get("outTree");
  TTree * exclChTree = (TTree*) mcFile->Get("exclChTree");

 
  // Tell the reader which branches are interesting
  mcTree->SetBranchStatus("*", 0);
  mcTree->SetBranchStatus("Ener2", 1);
  mcTree->SetBranchStatus("ESum2", 1);
  mcTree->SetBranchStatus("Multiplicity", 1);
  mcTree->SetBranchStatus("MultipletIndex", 1);
  mcTree->SetBranchStatus("Layer", 1);
  // Get the branch addresses
  double Ener2_data, ESum2_data;
  Short_t Multiplicity_data;
  Short_t MultipletIndex_data;
  int Channel_data, Detector_data, Layer_data;
  
  mcTree->SetBranchAddress("Ener2", &Ener2_data);
  mcTree->SetBranchAddress("ESum2", &ESum2_data);
  mcTree->SetBranchAddress("Multiplicity", &Multiplicity_data);
  mcTree->SetBranchAddress("MultipletIndex", &MultipletIndex_data);
  mcTree->SetBranchAddress("Layer", &Layer_data);
  int Dataset_data;
  bool Included_data;
  exclChTree->SetBranchAddress("Dataset", &Dataset_data);
  exclChTree->SetBranchAddress("Included", &Included_data);


  TTree * outTree = new TTree("outTree", "outTree");
  TTree * outexclChTree = new TTree("exclChTree", "exclChTree");
  TTree * outFriendTree = new TTree("shiftedEnergyTree", "shiftedEnergyTree");

  // Create new output tree
  double Ener2, ESum2;
  Short_t Multiplicity;
  Short_t MultipletIndex;
  int Channel, Detector, Layer;
  
  outTree->Branch("Ener2", &Ener2, "Ener2/D");
  outTree->Branch("ESum2", &ESum2, "ESum2/D");
  outTree->Branch("Multiplicity", &Multiplicity, "Multiplicity/S");
  outTree->Branch("MultipletIndex", &MultipletIndex, "MultipletIndex/S");
  outTree->Branch("Layer", &Layer, "Layer/I");

  // create output exclChTree
  int Dataset;
  bool Included;
  outexclChTree->Branch("Dataset", &Dataset, "Dataset/I");
  outexclChTree->Branch("Included", &Included, "Included/O");
  
  // Create new friend tree for the new output MC files
  double Ener2_shift_up, Ener2_shift_down;
  double ESum2_shift_up = 0;
  double ESum2_shift_down = 0 ;
  
  // Make the new friend TTree
  outFriendTree->Branch("Ener2_shift_up", &Ener2_shift_up, "Ener2_shift_up/D");
  outFriendTree->Branch("Ener2_shift_down", &Ener2_shift_down, "Ener2_shift_down/D");
  outFriendTree->Branch("ESum2_shift_up", &ESum2_shift_up, "ESum2_shift_up/D");
  outFriendTree->Branch("ESum2_shift_down", &ESum2_shift_down, "ESum2_shift_down/D");
  
  // Store the energies in a vector so that we can deal with Multiplicity > 1 events
  vector<double> energy_shift_up, energy_shift_down;
  
  bool previously_edited = false;

  // Get the starting index
  Long64_t total_entries = mcTree->GetEntries();
  Long64_t start_entry = std::floor(((double)(split-1) / number_of_splits) * total_entries);
  Long64_t end_entry = std::floor(((double)split / number_of_splits) * total_entries) - 1;

  // Find a good starting entry after the starting index (multiplicity == 1)
  bool found_start = false;
  bool found_end = false;
  Long64_t e = start_entry;
  while (!found_start)
    {
      mcTree->GetEntry(e);
      if (int(Multiplicity_data) == 1) {
	found_start = true;
	start_entry = e;
      }
      else e++;
    }
  // Find a good ending entry after the starting index (multiplicity == 1)
  e = end_entry;
  while (!found_end)
    {
      mcTree->GetEntry(e);
      if (int(Multiplicity_data) == 1) {
	found_end = true;
	end_entry = e;
      }
      else e--;
    }
  
  // Grab the events in the file
  e = start_entry;
  double percent;
  int i;
  while (e <= end_entry) {

    // Tell the user how far we are through the MC file every 10k events
    percent = (double(e - start_entry) / (end_entry - start_entry)) * 100.0;
    if ((e - start_entry) % 10000 == 0 | e == end_entry) {
      std::cout << "Entries read: " << (e - start_entry) << " of " << (end_entry - start_entry) << " | " << std::fixed <<  std::setprecision(4) << percent << "% complete" << "\r" << std::flush;
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
	  if (i + 1 >= int(Multiplicity_data))
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
    
    // Fill the outtree for each
    Ener2 = Ener2_data;
    ESum2 = ESum2_data;
    Multiplicity = Multiplicity_data;
    MultipletIndex = MultipletIndex_data;
Layer = Layer_data;
    Dataset = Dataset_data;
    Included_data = Included_data;

    outTree->Fill();
    outexclChTree->Fill();
    outFriendTree->Fill();
    // Go to next event
    e++;
  }
  std::cout << "\n";
  mcFile->Close();
  // Write the new output MC file
  outTree->Write();
  outexclChTree->Write();
  outFriendTree->Write();
  outFile->Close();
  delete outFile;
  delete mcFile;
}


int main(int argc, char **argv) {
  
  if (argc !=5) {
    std::cout << "usage: ./AddCalibrationFunction MC_file.root split number_of_splits directory" << std::endl;
    exit(1);
  }
  
  std::string mc_file = argv[1]; // The mc file is passed through the arguments of the function by the python script
  stringstream ss;
  ss << argv[2];
  int split;
  ss >> split; // the split that is considered here
  ss.str("");
  ss.clear();
  ss << argv[3];
  int number_of_splits;
  ss >> number_of_splits;  // the total number of splits
  ss.str("");
  ss.clear();
  std::string directory = argv[4];
  // The file containing the energy shifts
  std::string shifted_filename = "EnergyShift.txt";

  // Read the file and record the data as a map
  // Format looks like
  // || energy | shifted energy ||
  std::pair<std::map<double, double>, std::map<double, double> >  energy_shift_map_pair = ReadShift(shifted_filename);
  std::map<double, double> energy_shift_up = energy_shift_map_pair.first;
  std::map<double, double> energy_shift_down = energy_shift_map_pair.second;
  
  AddCalibrationFunction(mc_file, directory, energy_shift_up, energy_shift_down, split, number_of_splits);
  
}
