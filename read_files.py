# run with python 3

import argparse
import os
import sys

# read arguments
def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mc_list", help="the path to the MC_List.txt file")
    parser.add_argument("-i", "--iteration", type=int, help="the iteration of the file to run as a job")
    args = parser.parse_args()

    return args
    
# read in file to run on
def read_filename_and_path(mc_list, line_number):
    mc_file = open(mc_list, "r")
    i = 0
    for line in mc_file:
        #print(line)
        i = i + 1
        if (i == 1):
            directory = line.rstrip("\n")
        if (i == line_number):
            file_to_read = line.rstrip("\n")
    if i < line_number:
        print("line requested longer than file")
        sys.exit(1)
    path_to_read = directory + "/" + file_to_read
    return path_to_read

args = read_args()

mc_file = read_filename_and_path(args.mc_list, args.iteration)

os.system("./AddCalibrationFunction " + mc_file)
