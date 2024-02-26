import os
import h5py
import subprocess
import sys
# def check_subdirectories_for_file(root_dir, target_file):
#    directories_without_file = []
#    for dirpath, _, _ in os.walk(root_dir):
#        if os.path.isdir(dirpath):
#            files_in_dir = os.listdir(dirpath)
#            if target_file not in files_in_dir:
#                directories_without_file.append(dirpath)
#    return directories_without_file
#
# root_directory = './DFN_result'
# file_to_check = 'log.txt'
#
# directories_without_file = check_subdirectories_for_file(root_directory, file_to_check)
#
# if not directories_without_file:
#    print("All subdirectories contain '{}' file.".format(file_to_check))
# else:
#    print("The following directories do not contain '{}' file:".format(file_to_check))
#    for directory in directories_without_file:
#        print(directory)


def get_subdirectories(path):
    subdirectories = []
    for item in os.listdir(path):
        item_path = os.path.join(path, item)
        if os.path.isdir(item_path):
            subdirectories.append(item)
    return subdirectories


def file_exists(file_path):
    return os.path.exists(file_path)


def dataset_exists(file_path, dataset_name):
    exists = False
    with h5py.File(file_path, 'r') as file:
        if dataset_name in file:
            exists = True
    return exists


def call_bash_command(command):
    subprocess.run(command, shell=True)


def show_log(pathname):  # directory_path+"/"+subdir
    call_bash_command("cat " + pathname)
    ifcontinue = input("If continue? (input 0 or 1): ")
    if (ifcontinue == '0'):
        sys.exit()

# Replace 'path_to_your_directory' with the path to the directory you want to search
directory_path = './DFN_result_L_0030'

subdirectories = get_subdirectories(directory_path)

print("there are", len(subdirectories), "directories in", directory_path)

NumUnFinished = 0
for subdir in subdirectories:
    dirhdf5 = directory_path+"/"+subdir + "/Data.h5"
    if (not file_exists(dirhdf5)):
        print(dirhdf5, "does not exist")
        #show_log(directory_path+"/"+subdir + "/log.txt")
        NumUnFinished = NumUnFinished + 1
        continue
    exist_y = False
    if dataset_exists(dirhdf5, "Conn"):
        if dataset_exists(dirhdf5, "Permeab"):
            if dataset_exists(dirhdf5, "NumFiniteElements"):
                if dataset_exists(dirhdf5, "DFNTime"):
                    if dataset_exists(dirhdf5, "MeshTime"):
                        if dataset_exists(dirhdf5, "FlowTime"):
                            if dataset_exists(dirhdf5, "NumFractures"):
                                exist_y = True
    if (not exist_y):
        print(directory_path+"/"+subdir, "has not been finished")
        NumUnFinished = NumUnFinished + 1
        #show_log(directory_path+"/"+subdir + "/log.txt")
print("there are", NumUnFinished, "unfinished tasks")