import os
import argparse

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


def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mc_list", help="the path to the MC_List.txt file")
    args = parser.parse_args()

    return args

args = read_args()

mc_file_list = read_filename_and_path(args.mc_list)

i = 0
for mc_file in mc_file_list:

    print ("hadd file %s of %s" %(i, len(mc_file_list)))
    directory = mc_file.replace(".root", "")

    command = "hadd " + directory + "_adjusted.root " + directory + "/*"
    rm_command = "rm " + directory + "/*"


    print (command)
    os.system(command)
    #print (rm_command)

    
