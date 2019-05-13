# 2nu_CalibrationSystematic #
Scripts to add in a calibration systematic to the MC for the 2nu JAGS fit.
Note that the scripts are somewhat disk-intensive as it requires the MC files to be doubled in size, although some added cleverness can reduce this (and certainly post-operation deletions can also suffice).

## Steps ##

* Get the systematic information in a file named `EnergyShift.txt` with data organized into tab-delimited `Energy|shift|shift 1-sigma`
* Create a list of MC files to be edited (first line directory path, rest of lines mc filenames)
* Run `shift_single_file.pbs` to add two trees to each of the MC files
  * This is done in parts to speed up the computing time. As a result, there are some caveats and introduced biases in how it's done, but a reasonable splitting <100 should be totally fine.
  * Make sure to keep the number of jobs the same as the number of splittings, namely, `read_single_file_split.py -s 10` and `#PBS -t 1-10%10` should match
* Run `hadd.py` to compile all the output files into a new `_adjusted.root` file
* When running the JAGS scripts, edit `Parameters*.txt` to take in the new trees for the energy variables in the MC
  * `# Tree names\n` followed by `TreeName_MC outTree_down` or `TreeName_MC outTree_up`

