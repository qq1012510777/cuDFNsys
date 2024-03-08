import os


def check_file_existence(file_path):
    #print("--", file_path, ", ", os.path.exists(file_path))
    return os.path.exists(file_path)


def check_keyword_in_file(file_path, keyword):
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if keyword in line:
                    return True
        return False
    except FileNotFoundError:
        return False


current_directory = os.getcwd()

NUM_MC = 400

for i in range(1, NUM_MC+1):
    file_path = current_directory + "/MC_" + str(i).zfill(6)
    if ((not check_file_existence(file_path + "/Class_DFN.h5")) or (not check_file_existence(file_path + "/Class_MESH.h5"))
            or (not check_file_existence(file_path + "/Class_FLOW.h5"))):
        print("error in MC_" + str(i).zfill(6))

for i in range(1, NUM_MC+1):
    if (check_keyword_in_file(current_directory + "/log_MC_" + str(i+1).zfill(6), "Warning")):
        print("log error in MC_" + str(i).zfill(6))
