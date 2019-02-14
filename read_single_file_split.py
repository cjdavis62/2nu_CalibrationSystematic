# run with python 3

import argparse
import os
import sys

# read arguments
def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mc_list", help="the path to the MC_List.txt file")
    parser.add_argument("-s", "--split", type = int, help = "the number of ways to split the job")
    parser.add_argument("-n", "--number", type=int, help = "the number of this splitting")
    args = parser.parse_args()

    return args
    
# read in file to run on
def read_filename_and_path(mc_list):
    mc_file = open(mc_list, "r")
    file_list = []
    i = 0
    for line in mc_file:
        i = i + 1
        if (i == 1):
            directory = line.rstrip("\n")
        if (i != 1):
            file_to_read = line.rstrip("\n")
            path_to_read = directory + "/" + file_to_read
            file_list.append(path_to_read)
    return file_list

args = read_args()

mc_file_list = read_filename_and_path(args.mc_list)

for mc_file in mc_file_list:
    directory = mc_file.replace(".root","")
    print(mc_file)
    print (directory)
    try:
        os.mkdir(directory)
    except FileExistsError:
        pass
    os.system("./AddCalibrationFunction_doubletree " + mc_file + " " + str(args.number) + " " + str(args.split) + " " + directory)
