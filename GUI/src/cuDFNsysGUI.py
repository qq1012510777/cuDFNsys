from tkinter import *
import math
import subprocess
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import sys
import shutil


def generate_or_empty_file(file_name):
    # Open the file in write mode ('w'), which truncates the file if it exists or creates a new empty file if it doesn't
    with open(file_name, "w"):
        pass  # Using pass since we don't need to perform any write operation


def add_line_to_file(file_name, line):
    with open(file_name, "a") as file:  # Open the file in append mode ('a')
        file.write(line + "\n")


current_directory = os.getcwd()
# current_directory = os.path.dirname(os.path.abspath(__file__))
# current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
print("Current working directory:", current_directory)

cuDFNsys_root_dir = ""
cuDFNsys_GUI_path = ""


def DFN_input_para():
    # ----------------use existing csv
    def UseExistingCSV():
        def GenDFN_here():
            csv_name = entry1.get()
            seed_random = entry2.get()
            # print("CSV Name:", csv_name)
            try:
                # Replace "your_cpp_executable" with the path to your C++ executable
                cpp_executable = (
                    cuDFNsys_GUI_path + "/DFN_Gen " + csv_name + " " + seed_random
                )

                process = subprocess.Popen(cpp_executable, stdout=subprocess.PIPE,
                                           stderr=subprocess.STDOUT, universal_newlines=True, shell=True)
                # Read and print the output line by line
                for line in process.stdout:
                    print(line, end="")
                # Wait for the process to finish
                return_code = process.wait()

                new_window.destroy()
                if (return_code != 0):
                    raise subprocess.CalledProcessError(
                        return_code, cpp_executable, "Error")
                # ----------
                success_window = Toplevel()
                success_window.title("DFN generated")
                success_label = Label(success_window, text="Generated a DFN!")
                success_label.pack()
            except subprocess.CalledProcessError as e:
                print("Error running DFN_Gen:", e)
                new_window.destroy()
                # ----------
                success_window = Toplevel()
                success_window.title("Error")
                success_label = Label(success_window, text="Error happens!")
                success_label.pack()

        new_window = Toplevel()
        new_window.title("Input csv name")
        new_window.geometry("300x200")

        text = "Please input the name of the .csv file (without `.csv`):"
        label1 = Label(new_window, text=text, wraplength=250)
        label1.pack()
        entry1 = Entry(new_window)

        entry1.pack()

        text_ii = "Input the random seed (leave this blank for using the current time as a seed):"
        label2 = Label(new_window, text=text_ii, wraplength=250)
        label2.pack()
        entry2 = Entry(new_window)

        button = Button(new_window, text="Generate DFN", command=(GenDFN_here))

        entry2.pack()
        button.pack()

    # --------------------- determinsitic csv ----------------
    def GenDeterministicFractures():
        global Fracture_count
        global NumFracures_global
        Fracture_count = 1
        NumFracures_global = 1

        def GetDomainPara():
            global Fracture_count
            global NumFracures_global

            def writeFractures():
                global Fracture_count
                global NumFracures_global
                NormVec = entry_s1.get()
                Center_c = entry_s2.get()
                Radius_c = entry_s3.get()
                Beta_c = entry_s4.get()
                Gamma_c = entry_s5.get()
                AVertex = entry_s6.get()
                if not AVertex:
                    AVertex = ",,,"
                file_name = current_directory + "/DeterministicFracs"
                # print("Fracture_count: ", Fracture_count)
                add_line_to_file(
                    file_name + ".csv",
                    "Fracture_"
                    + str(Fracture_count)
                    + ", "
                    + NormVec
                    + ","
                    + Center_c
                    + ","
                    + Radius_c
                    + ","
                    + Beta_c
                    + ","
                    + Gamma_c
                    + ","
                    + AVertex
                    + ",",
                )
                button_s1.pack_forget()
                text_button_s2 = "Continue"
                if Fracture_count == NumFracures_global:
                    text_button_s2 = "Generate DFN"
                button_s2 = Button(
                    add_frac_window[Fracture_count - 1],
                    text=text_button_s2,
                    command=(GetDomainPara),
                )
                button_s2.pack()
                Fracture_count = Fracture_count + 1

                # ---------------------

            if Fracture_count == 1:
                DomainSizeX = entry1.get()
                DomainDimensionRatio = entry2.get()
                if not DomainDimensionRatio:
                    DomainDimensionRatio = "1, 1, 1"
                Percolation_direction = entry3.get()
                if not Percolation_direction:
                    Percolation_direction = "2"
                NumFractures = entry4.get()
                NumFracures_global = int(NumFractures)
                file_name = current_directory + "/DeterministicFracs"
                add_line_to_file(file_name + ".csv",
                                 "DomainSizeX," + DomainSizeX + ",")
                add_line_to_file(
                    file_name + ".csv",
                    "DomainDimensionRatio," + DomainDimensionRatio + ",",
                )
                add_line_to_file(
                    file_name + ".csv",
                    "Percolation_direction," + Percolation_direction + ",",
                )
                add_line_to_file(
                    file_name + ".csv", "NumFractures," + NumFractures + ","
                )
                new_window.destroy()
            # ----- generate a new window to add frac
            # print("GetDomainPara: Fracture_count: ",
            #       Fracture_count, NumFracures_global)
            if Fracture_count > 1:
                add_frac_window[Fracture_count - 2].destroy()
            if Fracture_count > NumFracures_global:
                try:
                    # Replace "your_cpp_executable" with the path to your C++ executable
                    file_name = current_directory + "/DeterministicFracs"
                    cpp_executable = cuDFNsys_GUI_path + "/DFN_Gen " + file_name
                    # output = subprocess.check_output(
                    #    [cpp_executable], shell=True)
                    # print("DFN_Gen output:")
                    # print(output.decode("utf-8"))  # Decode bytes to string
                    process = subprocess.Popen(
                        cpp_executable, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, shell=True)
                    # Read and print the output line by line
                    for line in process.stdout:
                        print(line, end="")
                    # Wait for the process to finish
                    return_code = process.wait()
                    new_window.destroy()
                    if (return_code != 0):
                        raise subprocess.CalledProcessError(
                            return_code, cpp_executable, "Error")
                    # ----------
                    success_window = Toplevel()
                    success_window.title("DFN generated")
                    success_label = Label(
                        success_window, text="Generated a DFN!")
                    success_label.pack()
                except subprocess.CalledProcessError as e:
                    print("Error running DFN_Gen:", e)
                    new_window.destroy()
                    # ----------
                    success_window = Toplevel()
                    success_window.title("Error")
                    success_label = Label(
                        success_window, text="Error happens!")
                    success_label.pack()
                return
            New_ADDFRAC_window = Toplevel()
            New_ADDFRAC_window.title("Add fracture No. " + str(Fracture_count))
            New_ADDFRAC_window.geometry("300x390")

            text_s1 = "Normal vector (normalized) (e.g., 0, 0, 1):"
            label_s1 = Label(New_ADDFRAC_window,
                             text=text_s1, wraplength=425000)
            label_s1.pack()
            entry_s1 = Entry(New_ADDFRAC_window)
            entry_s1.pack()

            text_s2 = "Fracture center (e.g., 0, 0, 0):"
            label_s2 = Label(New_ADDFRAC_window, text=text_s2, wraplength=250)
            label_s2.pack()
            entry_s2 = Entry(New_ADDFRAC_window)
            entry_s2.pack()

            text_s3 = "Fracture radius R (e.g., 40):"
            label_s3 = Label(New_ADDFRAC_window, text=text_s3, wraplength=250)
            label_s3.pack()
            entry_s3 = Entry(New_ADDFRAC_window)
            entry_s3.pack()

            text_s4 = (
                r"Beta (note that k_f = b_f^3 / 12 =  [(Gamma * R ^ Beta)] ^ 3 / 12):"
            )
            label_s4 = Label(New_ADDFRAC_window, text=text_s4, wraplength=250)
            label_s4.pack()
            entry_s4 = Entry(New_ADDFRAC_window)
            entry_s4.pack()

            text_s5 = "Gamma:"
            label_s5 = Label(New_ADDFRAC_window, text=text_s5, wraplength=250)
            label_s5.pack()
            entry_s5 = Entry(New_ADDFRAC_window)
            entry_s5.pack()

            text_s6 = "A vertex of the square fracture, e.g., 1, 1, 1 (leave it blank for random selections):"
            label_s6 = Label(New_ADDFRAC_window, text=text_s6, wraplength=250)
            label_s6.pack()
            entry_s6 = Entry(New_ADDFRAC_window)
            entry_s6.pack()

            text_button_s1 = "Add fracture No. " + str(Fracture_count)
            button_s1 = Button(
                New_ADDFRAC_window, text=text_button_s1, command=(
                    writeFractures)
            )
            button_s1.pack()
            add_frac_window.append(New_ADDFRAC_window)

        new_window = Toplevel()
        new_window.title("Generate deterministic fractures")
        new_window.geometry("300x300")

        file_name = current_directory + "/DeterministicFracs"
        generate_or_empty_file(file_name + ".csv")
        add_line_to_file(file_name + ".csv", "IfStochastic,0,")

        text1 = "Domain size in the x direction (e.g., 30):"
        label1 = Label(new_window, text=text1, wraplength=250)
        label1.pack()
        entry1 = Entry(new_window)
        entry1.pack()

        text2 = "Domain dimension ratio (e.g., 1, 2, 3)(leave it blank for 1, 1, 1):"
        label2 = Label(new_window, text=text2, wraplength=250)
        label2.pack()
        entry2 = Entry(new_window)
        entry2.pack()

        text3 = "Percolation direction (input: 0 for x, or 1 for y, or 2 for z)(leave it blank for 2):"
        label3 = Label(new_window, text=text3, wraplength=250)
        label3.pack()
        entry3 = Entry(new_window)
        entry3.pack()

        text4 = "Number of fractures (e.g., 4):"
        label4 = Label(new_window, text=text4, wraplength=250)
        label4.pack()
        entry4 = Entry(new_window)
        entry4.pack()

        add_frac_window = []

        button1 = Button(new_window, text="Continue", command=(GetDomainPara))
        button1.pack()

    # ---------------------stochastic csv
    def GenStochasticFractures():
        generate_or_empty_file(current_directory + "/StochasticFracs.csv")
        add_line_to_file(current_directory +
                         "/StochasticFracs.csv", "IfStochastic,1,")
        NumFractureGroup = 1

        def RecordFirstFourData():
            global NumFractureGroup
            DomainSizeX = entry1.get()
            DomainDimensionRatio = entry2.get()
            if not DomainDimensionRatio:
                DomainDimensionRatio = "1,1,1"
            Percolation_direction = entry3.get()
            if not Percolation_direction:
                Percolation_direction = "2"
            NumFractureGroup = entry4.get()
            if not NumFractureGroup:
                NumFractureGroup = "1"
            add_line_to_file(
                current_directory + "/StochasticFracs.csv",
                "DomainSizeX," + str(DomainSizeX) + ",",
            )
            add_line_to_file(
                current_directory + "/StochasticFracs.csv",
                "DomainDimensionRatio," + str(DomainDimensionRatio) + ",",
            )
            add_line_to_file(
                current_directory + "/StochasticFracs.csv",
                "Percolation_direction," + str(Percolation_direction) + ",",
            )
            add_line_to_file(
                current_directory + "/StochasticFracs.csv",
                "NumFractureGroups," + str(NumFractureGroup) + ",",
            )
            for widget in StochasticWindow_1.winfo_children():
                widget.destroy()
            # ----------------------------------
            global Group_present
            Group_present = 1
            global NumFracs
            global Kappa
            global MeanOri
            global ModeSize
            global ParaSize
            global Beta
            global Gamma
            NumFracs = ""
            Kappa = ""
            MeanOri = ""
            ModeSize = ""
            ParaSize = ""
            Beta = ""
            Gamma = ""

            def DataOfEachGroup():
                global Group_present
                global NumFractureGroup
                global NumFracs
                global Kappa
                global MeanOri
                global ModeSize
                global ParaSize
                global Beta
                global Gamma
                global entry_k1
                global entry_k2
                global entry_k3
                global entry_k4
                global entry_k5
                global entry_k6
                global entry_k7
                if Group_present > 1:
                    NumFracs_tt = entry_k1.get()
                    Kappa_tt = entry_k2.get()
                    if not Kappa_tt:
                        Kappa_tt = "0"
                    MeanOri_tt = entry_k3.get()
                    if not MeanOri_tt:
                        MeanOri_tt = "0, 0, 1"
                    ModeSize_tt = entry_k4.get()
                    ParaSize_tt = entry_k5.get()
                    Beta_tt = entry_k6.get()
                    Gamma_tt = entry_k7.get()

                    NumFracs = NumFracs + (NumFracs_tt) + ","
                    Kappa = Kappa + (Kappa_tt) + ","
                    MeanOri = MeanOri + (MeanOri_tt) + ","
                    ModeSize = ModeSize + (ModeSize_tt) + ","
                    ParaSize = ParaSize + (ParaSize_tt) + ","
                    Beta = Beta + (Beta_tt) + ","
                    Gamma = Gamma + (Gamma_tt) + ","
                    for widget in StochasticWindow_1.winfo_children():
                        widget.destroy()

                if Group_present > int(NumFractureGroup):

                    add_line_to_file(
                        current_directory + "/StochasticFracs.csv",
                        "NumFractureEachGroup," + NumFracs,
                    )

                    add_line_to_file(
                        current_directory + "/StochasticFracs.csv",
                        "KappaValues," + Kappa,
                    )

                    add_line_to_file(
                        current_directory + "/StochasticFracs.csv",
                        "MeanOrientationOfFisherDistribution," + MeanOri,
                    )

                    add_line_to_file(
                        current_directory + "/StochasticFracs.csv",
                        "ModeOfSizeDistribution," + ModeSize,
                    )

                    add_line_to_file(
                        current_directory + "/StochasticFracs.csv",
                        "SizeDistributionParameters," + ParaSize,
                    )

                    add_line_to_file(
                        current_directory + "/StochasticFracs.csv",
                        "Beta," + Beta,
                    )

                    add_line_to_file(
                        current_directory + "/StochasticFracs.csv",
                        "Gamma," + Gamma,
                    )

                    StochasticWindow_1.title("Random seed:")

                    text_k1 = (
                        "Random seeds (leave it blank for using the current as a seed):"
                    )
                    label_k1 = Label(StochasticWindow_1,
                                     text=text_k1, wraplength=250)
                    label_k1.pack()
                    entry_k1 = Entry(StochasticWindow_1)
                    entry_k1.pack()

                    global randomseed_sd
                    randomseed_sd = ""

                    def GetRandomseed():
                        global randomseed_sd
                        randomseed_sd = entry_k1.get()
                        try:
                            # Replace "your_cpp_executable" with the path to your C++ executable
                            cpp_executable = (
                                cuDFNsys_GUI_path
                                + "/DFN_Gen "
                                + current_directory
                                + "/StochasticFracs "
                                + randomseed_sd
                            )
                            # output = subprocess.check_output(
                            #    [cpp_executable], shell=True)
                            # print("DFN_Gen output:")
                            # print(output.decode("utf-8"))  # Decode bytes to string
                            process = subprocess.Popen(
                                cpp_executable, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, shell=True)
                            # Read and print the output line by line
                            for line in process.stdout:
                                print(line, end="")
                            # Wait for the process to finish
                            return_code = process.wait()
                            StochasticWindow_1.destroy()
                            if (return_code != 0):
                                raise subprocess.CalledProcessError(
                                    return_code, cpp_executable, "Error")
                            # ----------
                            success_window = Toplevel()
                            success_window.title("DFN generated")
                            success_label = Label(
                                success_window, text="Generated a DFN!")
                            success_label.pack()
                        except subprocess.CalledProcessError as e:
                            print("Error running DFN_Gen:", e)
                            StochasticWindow_1.destroy()
                            # ----------
                            success_window = Toplevel()
                            success_window.title("Error")
                            success_label = Label(
                                success_window, text="Error happens!")
                            success_label.pack()

                    button1 = Button(
                        StochasticWindow_1, text="Generate DFN", command=(GetRandomseed)
                    )
                    button1.pack()  

                    return

                # --------------
                StochasticWindow_1.title(
                    "Inputs of group No. " + str(Group_present))
                StochasticWindow_1.geometry("300x600")
                text_k1 = "Number of fractures:"
                label_k1 = Label(StochasticWindow_1,
                                 text=text_k1, wraplength=250)
                label_k1.pack()
                entry_k1 = Entry(StochasticWindow_1)
                entry_k1.pack()

                text_k2 = "Kappa value (leave it blank for 0):"
                label_k2 = Label(StochasticWindow_1,
                                 text=text_k2, wraplength=250)
                label_k2.pack()
                entry_k2 = Entry(StochasticWindow_1)
                entry_k2.pack()

                text_k3 = "Mean orientation (e.g., 0, 1, 0)(normalized)(leave it blank for 0, 0, 1):"
                label_k3 = Label(StochasticWindow_1,
                                 text=text_k3, wraplength=250)
                label_k3.pack()
                entry_k3 = Entry(StochasticWindow_1)
                entry_k3.pack()

                text_k4 = "Mode of size distribution (0: power-law, 1: lognormal, 2: uniform, 3: mono-sized):"
                label_k4 = Label(StochasticWindow_1,
                                 text=text_k4, wraplength=250)
                label_k4.pack()
                entry_k4 = Entry(StochasticWindow_1)
                entry_k4.pack()

                text_k5 = "Parameters of size distribution (note that there must be four values separatered by commas, e.g., if you input 2 in the above entry for uniform distribution, you can input 1, 10, 0, 0, meaning minimum and maximum are 1 and 10, more details are in Manual):"
                label_k5 = Label(StochasticWindow_1,
                                 text=text_k5, wraplength=250)
                label_k5.pack()
                entry_k5 = Entry(StochasticWindow_1)
                entry_k5.pack()

                text_k6 = "Beta (note that k_f = b_f^3 / 12 =  [(Gamma * R ^ Beta)] ^ 3 / 12):"
                label_k6 = Label(StochasticWindow_1,
                                 text=text_k6, wraplength=250)
                label_k6.pack()
                entry_k6 = Entry(StochasticWindow_1)
                entry_k6.pack()

                text_k7 = "Gamma:"
                label_k7 = Label(StochasticWindow_1,
                                 text=text_k7, wraplength=250)
                label_k7.pack()
                entry_k7 = Entry(StochasticWindow_1)
                entry_k7.pack()

                textButton = "Continue"
                button1 = Button(
                    StochasticWindow_1, text=textButton, command=(
                        DataOfEachGroup)
                )
                button1.pack()
                Group_present = Group_present + 1

            DataOfEachGroup()

        StochasticWindow_1 = Toplevel()
        StochasticWindow_1.title("Generate stochastic fractures")
        StochasticWindow_1.geometry("400x300")

        text1 = "Domain size in the x direction (e.g., 30):"
        label1 = Label(StochasticWindow_1, text=text1, wraplength=250)
        label1.pack()
        entry1 = Entry(StochasticWindow_1)
        entry1.pack()

        text2 = "Domain dimension ratio (e.g., 1, 2, 3)(leave it blank for 1, 1, 1):"
        label2 = Label(StochasticWindow_1, text=text2, wraplength=250)
        label2.pack()
        entry2 = Entry(StochasticWindow_1)
        entry2.pack()

        text3 = "Percolation direction (input: 0 for x, or 1 for y, or 2 for z)(leave it blank for 2):"
        label3 = Label(StochasticWindow_1, text=text3, wraplength=250)
        label3.pack()
        entry3 = Entry(StochasticWindow_1)
        entry3.pack()

        text4 = "Number of fracture groups (e.g., 2)(leave it blank for 1 group):"
        label4 = Label(StochasticWindow_1, text=text4, wraplength=250)
        label4.pack()
        entry4 = Entry(StochasticWindow_1)
        entry4.pack()

        button1 = Button(
            StochasticWindow_1, text="Continue", command=(RecordFirstFourData)
        )
        button1.pack()
        # print(NumFractureGroup)

    def VisualizeDFN():
        def UseMatPlotLib():
            f = h5py.File(current_directory + "/Class_DFN.h5")
            NumFracs = np.array(f["NumFractures"][0])
            Lm = np.array(f["L"][0])
            DomainDimensionRatio = np.array(f["DomainDimensionRatio"][:])
            # print(NumFracs)
            Verts = np.zeros((NumFracs, 4, 3))
            # Faces = []

            # Create a 3D plot
            fig = plt.figure()
            ax = fig.add_subplot(111, projection="3d")
            MaxR = 0

            for i in range(0, NumFracs):
                vsd = np.array(f["Fracture_" + str(i + 1) + "/Verts3D"][:])
                vsd = np.transpose(vsd)
                Verts[i, :, :] = vsd
                R_rr = np.linalg.norm(vsd[1, :] - vsd[3, :])
                if R_rr > MaxR:
                    MaxR = R_rr
                # Faces.append([0, 1, 2, 3])
            for verts in Verts:
                poly = Poly3DCollection(
                    [verts], facecolors="red", linewidths=1, edgecolors="k", alpha=1
                )
                ax.add_collection3d(poly)

            # Define coordinates of cuboid vertices
            vertices = np.array(
                [
                    [-DomainDimensionRatio[0] * Lm / 2, -DomainDimensionRatio[1]
                        * Lm / 2, -DomainDimensionRatio[2] * Lm / 2],
                    [-DomainDimensionRatio[0] * Lm / 2, -DomainDimensionRatio[1]
                        * Lm / 2, DomainDimensionRatio[2] * Lm / 2],
                    [-DomainDimensionRatio[0] * Lm / 2, DomainDimensionRatio[1]
                        * Lm / 2, -DomainDimensionRatio[2] * Lm / 2],
                    [-DomainDimensionRatio[0] * Lm / 2, DomainDimensionRatio[1]
                        * Lm / 2, DomainDimensionRatio[2] * Lm / 2],
                    [DomainDimensionRatio[0] * Lm / 2, -DomainDimensionRatio[1]
                        * Lm / 2, -DomainDimensionRatio[2] * Lm / 2],
                    [DomainDimensionRatio[0] * Lm / 2, -DomainDimensionRatio[1]
                        * Lm / 2, DomainDimensionRatio[2] * Lm / 2],
                    [DomainDimensionRatio[0] * Lm / 2, DomainDimensionRatio[1]
                        * Lm / 2, -DomainDimensionRatio[2] * Lm / 2],
                    [DomainDimensionRatio[0] * Lm / 2, DomainDimensionRatio[1]
                        * Lm / 2, DomainDimensionRatio[2] * Lm / 2],
                ]
            )
            # Define edges of cuboid
            edges = [
                [vertices[0], vertices[1], vertices[3], vertices[2], vertices[0]],
                [vertices[4], vertices[5], vertices[7], vertices[6], vertices[4]],
                [vertices[0], vertices[4]],
                [vertices[1], vertices[5]],
                [vertices[2], vertices[6]],
                [vertices[3], vertices[7]],
            ]
            for edge in edges:
                ax.add_collection3d(Line3DCollection(
                    [edge], color="black", linewidths=1))

            ax.set_xlim(
                [
                    -DomainDimensionRatio[0] * Lm / 2.0 - R_rr / 2.0,
                    DomainDimensionRatio[0] * Lm / 2.0 + R_rr / 2.0,
                ]
            )
            ax.set_ylim(
                [
                    -DomainDimensionRatio[1] * Lm / 2.0 - R_rr / 2.0,
                    DomainDimensionRatio[1] * Lm / 2.0 + R_rr / 2.0,
                ]
            )
            ax.set_zlim(
                [
                    -DomainDimensionRatio[2] * Lm / 2.0 - R_rr / 2.0,
                    DomainDimensionRatio[2] * Lm / 2.0 + R_rr / 2.0,
                ]
            )
            ax.set_box_aspect([1, 1, 1])
            ax.set_xlabel("X (m)")
            ax.set_ylabel("Y (m)")
            ax.set_zlabel("Z (m)")
            f.close()
            plt.show()
            visualizeDFN_window.destroy()

        def UseMayavi():
            try:
                # Replace "your_cpp_executable" with the path to your C++ executable
                cpp_executable = "python3 " + current_directory + "/DFN_VISUAL.py"
                # output = subprocess.check_output([cpp_executable], shell=True)
                # print("DFN_Gen output:")
                # print(output.decode("utf-8"))  # Decode bytes to string
                process = subprocess.Popen(cpp_executable, stdout=subprocess.PIPE,
                                           stderr=subprocess.STDOUT, universal_newlines=True, shell=True)
                # Read and print the output line by line
                for line in process.stdout:
                    print(line, end="")
                # Wait for the process to finish
                return_code = process.wait()

                visualizeDFN_window.destroy()
                # ----------
                if (return_code != 0):
                    raise subprocess.CalledProcessError(
                        return_code, cpp_executable, "Error")

            except subprocess.CalledProcessError as e:
                print("Error running DFN_Gen:", e)
                visualizeDFN_window.destroy()
                # ----------
            return

        visualizeDFN_window = Toplevel()
        visualizeDFN_window.title("Select DFN visualization methods:")
        visualizeDFN_window.geometry("300x70")

        btn1_t = Button(
            visualizeDFN_window, text="Use Matplotlib (not good)", command=UseMatPlotLib
        )
        btn1_t.pack(fill=X, expand=1)

        btn2_t = Button(
            visualizeDFN_window,
            text="Use Mayavi (require this library to be installed)",
            command=UseMayavi,
        )
        btn2_t.pack(fill=X, expand=1)

        return

    # -----------------DFN Gen window
    new_window = Toplevel()
    new_window.title("DFN generation inputs")
    new_window.geometry("300x130")

    btn1 = Button(
        new_window, text="Deterministic fractures", command=GenDeterministicFractures
    )
    btn1.pack(fill=X, expand=1)

    btn2 = Button(
        new_window, text="Stochastic fractures", command=GenStochasticFractures
    )
    btn2.pack(fill=X, expand=1)

    btn3 = Button(new_window, text="Use existing .csv file",
                  command=UseExistingCSV)
    btn3.pack(fill=X, expand=1)

    btn3 = Button(new_window, text="Visualize an existing DFN",
                  command=VisualizeDFN)
    btn3.pack(fill=X, expand=1)


def MeshOption():
    def DFN_mesh():
        def Gen_mesh():

            mixS = entry1.get()
            maxS = entry2.get()
            generate_or_empty_file(current_directory + "/MeshPara.csv")
            add_line_to_file(current_directory + "/MeshPara.csv",
                             "ExpectedMinimimGridSize, " + str(mixS) + ",")
            add_line_to_file(current_directory + "/MeshPara.csv",
                             "ExpectedMaximimGridSize, " + str(maxS) + ",")

            try:
                cpp_executable = (
                    cuDFNsys_GUI_path + "/DFN_Mesh " + current_directory + "/MeshPara"
                )
                # output = subprocess.check_output([cpp_executable], shell=True)
                # print("DFN_Mesh output:")
                # print(output.decode("utf-8"))  # Decode bytes to string

                process = subprocess.Popen(cpp_executable, stdout=subprocess.PIPE,
                                           stderr=subprocess.STDOUT, universal_newlines=True, shell=True)
                # Read and print the output line by line
                for line in process.stdout:
                    print(line, end="")

                # Wait for the process to finish
                return_code = process.wait()

                new_window.destroy()
                if (return_code != 0):
                    raise subprocess.CalledProcessError(
                        return_code, cpp_executable, "Error")
                # ----------
                success_window = Toplevel()
                success_window.title("DFN mesh generated")
                success_label = Label(
                    success_window, text="Generated a DFN mesh!")
                success_label.pack()
            except subprocess.CalledProcessError as e:
                print("Error running DFN_Gen:", e)
                new_window.destroy()
                # ----------
                success_window = Toplevel()
                success_window.title("Error")
                success_label = Label(success_window, text="Error happens!")
                success_label.pack()

        new_window = Toplevel()
        new_window.title("DFN mesh parameter:")
        new_window.geometry("300x130")

        text1 = "Expected minimum grid size:"
        label1 = Label(new_window, text=text1, wraplength=250)
        label1.pack()
        entry1 = Entry(new_window)
        entry1.pack()
        text2 = "Expected maximum grid size:"
        label2 = Label(new_window, text=text2, wraplength=250)
        label2.pack()
        entry2 = Entry(new_window)
        entry2.pack()

        btn1 = Button(
            new_window, text="Run mesh generation (take some time)", command=Gen_mesh)
        btn1.pack(fill=X, expand=1)

    def visualize_mesh():
        def matplotlibmesh():
            f = h5py.File(current_directory + "/Class_MESH.h5")
            Points = np.array(f["group_mesh/Coordinate3D"][:])
            Points = np.transpose(Points)

            Triangles = np.array(f["group_mesh/Element3D"][:])
            Triangles = np.transpose(Triangles)
            f.close()
            # print(Triangles)

            # Create a 3D plot
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            # Plot the mesh
            ax.plot_trisurf(Points[:, 0], Points[:, 1],
                            Points[:, 2], triangles=Triangles-1)

            f2 = h5py.File(current_directory + "/Class_DFN.h5")
            DomainDimensionRatio = np.array(f2["DomainDimensionRatio"][:])
            Lm = np.array(f2["L"][0])
            f2.close()

            ax.set_xlim(
                [
                    -DomainDimensionRatio[0] * Lm / 2.0,
                    DomainDimensionRatio[0] * Lm / 2.0,
                ]
            )
            ax.set_ylim(
                [
                    -DomainDimensionRatio[1] * Lm / 2.0,
                    DomainDimensionRatio[1] * Lm / 2.0,
                ]
            )
            ax.set_zlim(
                [
                    -DomainDimensionRatio[2] * Lm / 2.0,
                    DomainDimensionRatio[2] * Lm / 2.0,
                ]
            )
            vertices = np.array(
                [
                    [-DomainDimensionRatio[0] * Lm / 2, -DomainDimensionRatio[1]
                        * Lm / 2, -DomainDimensionRatio[2] * Lm / 2],
                    [-DomainDimensionRatio[0] * Lm / 2, -DomainDimensionRatio[1]
                        * Lm / 2, DomainDimensionRatio[2] * Lm / 2],
                    [-DomainDimensionRatio[0] * Lm / 2, DomainDimensionRatio[1]
                        * Lm / 2, -DomainDimensionRatio[2] * Lm / 2],
                    [-DomainDimensionRatio[0] * Lm / 2, DomainDimensionRatio[1]
                        * Lm / 2, DomainDimensionRatio[2] * Lm / 2],
                    [DomainDimensionRatio[0] * Lm / 2, -DomainDimensionRatio[1]
                        * Lm / 2, -DomainDimensionRatio[2] * Lm / 2],
                    [DomainDimensionRatio[0] * Lm / 2, -DomainDimensionRatio[1]
                        * Lm / 2, DomainDimensionRatio[2] * Lm / 2],
                    [DomainDimensionRatio[0] * Lm / 2, DomainDimensionRatio[1]
                        * Lm / 2, -DomainDimensionRatio[2] * Lm / 2],
                    [DomainDimensionRatio[0] * Lm / 2, DomainDimensionRatio[1]
                        * Lm / 2, DomainDimensionRatio[2] * Lm / 2],
                ]
            )
            # Define edges of cuboid
            edges = [
                [vertices[0], vertices[1], vertices[3], vertices[2], vertices[0]],
                [vertices[4], vertices[5], vertices[7], vertices[6], vertices[4]],
                [vertices[0], vertices[4]],
                [vertices[1], vertices[5]],
                [vertices[2], vertices[6]],
                [vertices[3], vertices[7]],
            ]
            for edge in edges:
                ax.add_collection3d(Line3DCollection(
                    [edge], color="black", linewidths=1))
            # Set labels
            ax.set_box_aspect([1, 1, 1])
            ax.set_xlabel('X (m)')
            ax.set_ylabel('Y (m)')
            ax.set_zlabel('Z (m)')

            plt.show()

        def mayavimesh():
            try:
                # Replace "your_cpp_executable" with the path to your C++ executable
                cpp_executable = "python3 " + current_directory + "/DFN_MESH_VISUAL.py"
                # output = subprocess.check_output([cpp_executable], shell=True)
                # print("DFN_Gen output:")
                # print(output.decode("utf-8"))  # Decode bytes to string
                # ----------
                process = subprocess.Popen(cpp_executable, stdout=subprocess.PIPE,
                                           stderr=subprocess.STDOUT, universal_newlines=True, shell=True)
                # Read and print the output line by line
                for line in process.stdout:
                    print(line, end="")
                # Wait for the process to finish
                return_code = process.wait()
                if (return_code != 0):
                    raise subprocess.CalledProcessError(
                        return_code, cpp_executable, "Error")
            except subprocess.CalledProcessError as e:
                print("Error running DFN_Gen:", e)
                # ----------

        for widget in mainWindowMesh.winfo_children():
            widget.destroy()
        mainWindowMesh.title("DFN mesh visualization:")
        btn1 = Button(mainWindowMesh, text="Use Matplotlib",
                      command=matplotlibmesh)
        btn1.pack(fill=X, expand=1)

        btn2 = Button(
            mainWindowMesh, text="Use Mayavi (require installed mayavi)", command=mayavimesh)
        btn2.pack(fill=X, expand=1)

    mainWindowMesh = Toplevel()
    mainWindowMesh.title("DFN mesh")
    mainWindowMesh.geometry("300x130")

    btn1 = Button(mainWindowMesh, text="Generate mesh", command=DFN_mesh)
    btn1.pack(fill=X, expand=1)

    btn2 = Button(mainWindowMesh, text="Visualize mesh",
                  command=visualize_mesh)
    btn2.pack(fill=X, expand=1)


def FlowOptionWindow():
    def FlowParaSet():
        def SloveFlow():
            generate_or_empty_file(current_directory + "/FlowPara.csv")

            add_line_to_file(current_directory + "/FlowPara.csv", "MuOverRhoG, " + str(entry1.get()) +
                             ",")
            add_line_to_file(current_directory + "/FlowPara.csv", "InletHead, " + str(entry2.get()) + ","
                             )
            add_line_to_file(current_directory + "/FlowPara.csv", "OutletHead, " + str(entry3.get()) +
                             ",")
            try:
                cpp_executable = (
                    cuDFNsys_GUI_path + "/DFN_Flow " + current_directory + "/FlowPara"
                )
                process = subprocess.Popen(cpp_executable, stdout=subprocess.PIPE,
                                           stderr=subprocess.STDOUT, universal_newlines=True, shell=True)
                for line in process.stdout:
                    print(line, end="")

                # Wait for the process to finish
                return_code = process.wait()
                FlowParaWindow.destroy()
                if (return_code != 0):
                    raise subprocess.CalledProcessError(
                        return_code, cpp_executable, "Error")
                # ----------
                success_window = Toplevel()
                success_window.title("DFN flow solver")
                success_label = Label(
                    success_window, text="Flow simulation finished")
                success_label.pack()
            except subprocess.CalledProcessError as e:
                print("Error running DFN_Flow:", e)
                FlowParaWindow.destroy()
                # ----------
                success_window = Toplevel()
                success_window.title("Error")
                success_label = Label(success_window, text="Error happens!")
                success_label.pack()
            return
        FlowParaWindow = Toplevel()
        FlowParaWindow.title("Input flow parameters")
        FlowParaWindow.geometry("400x300")

        text1 = r"Mu / (rho * g) (fluid viscosity over (density times gravity)):"
        label1 = Label(FlowParaWindow, text=text1, wraplength=250)
        label1.pack()
        entry1 = Entry(FlowParaWindow)
        entry1.pack()

        text2 = "Inlet head [L]:"
        label2 = Label(FlowParaWindow, text=text2, wraplength=250)
        label2.pack()
        entry2 = Entry(FlowParaWindow)
        entry2.pack()

        text3 = "Outlet head [L]:"
        label3 = Label(FlowParaWindow, text=text3, wraplength=250)
        label3.pack()
        entry3 = Entry(FlowParaWindow)
        entry3.pack()

        btn_flowpara = Button(
            FlowParaWindow, text="Solve flow", command=SloveFlow)
        btn_flowpara.pack()
        return

    def VisualizeFlowDFN():
        def MatplotlibVisualFlow():
            f = h5py.File(current_directory + "/Class_MESH.h5")
            Points = np.array(f["group_mesh/Coordinate3D"][:])
            Points = np.transpose(Points)

            Triangles = np.array(f["group_mesh/Element3D"][:])
            Triangles = np.transpose(Triangles)
            f.close()
            # print(Triangles)

            # Create a 3D plot
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            # Plot the mesh
            f3 = h5py.File(current_directory + "/Class_FLOW.h5")
            PressureElements = np.array(f3["PressureEles"][:])
            f3.close()
            surf = ax.plot_trisurf(Points[:, 0], Points[:, 1],
                                   Points[:, 2], triangles=Triangles-1, cmap='viridis', shade=False)
            # Set colors based on PressureElements values
            surf.set_array(PressureElements.ravel())
            fig.colorbar(surf, ax=ax, label='Hydraulic head')

            f2 = h5py.File(current_directory + "/Class_DFN.h5")
            DomainDimensionRatio = np.array(f2["DomainDimensionRatio"][:])
            Lm = np.array(f2["L"][0])
            f2.close()

            ax.set_xlim(
                [
                    -DomainDimensionRatio[0] * Lm / 2.0,
                    DomainDimensionRatio[0] * Lm / 2.0,
                ]
            )
            ax.set_ylim(
                [
                    -DomainDimensionRatio[1] * Lm / 2.0,
                    DomainDimensionRatio[1] * Lm / 2.0,
                ]
            )
            ax.set_zlim(
                [
                    -DomainDimensionRatio[2] * Lm / 2.0,
                    DomainDimensionRatio[2] * Lm / 2.0,
                ]
            )
            vertices = np.array(
                [
                    [-DomainDimensionRatio[0] * Lm / 2, -DomainDimensionRatio[1]
                        * Lm / 2, -DomainDimensionRatio[2] * Lm / 2],
                    [-DomainDimensionRatio[0] * Lm / 2, -DomainDimensionRatio[1]
                        * Lm / 2, DomainDimensionRatio[2] * Lm / 2],
                    [-DomainDimensionRatio[0] * Lm / 2, DomainDimensionRatio[1]
                        * Lm / 2, -DomainDimensionRatio[2] * Lm / 2],
                    [-DomainDimensionRatio[0] * Lm / 2, DomainDimensionRatio[1]
                        * Lm / 2, DomainDimensionRatio[2] * Lm / 2],
                    [DomainDimensionRatio[0] * Lm / 2, -DomainDimensionRatio[1]
                        * Lm / 2, -DomainDimensionRatio[2] * Lm / 2],
                    [DomainDimensionRatio[0] * Lm / 2, -DomainDimensionRatio[1]
                        * Lm / 2, DomainDimensionRatio[2] * Lm / 2],
                    [DomainDimensionRatio[0] * Lm / 2, DomainDimensionRatio[1]
                        * Lm / 2, -DomainDimensionRatio[2] * Lm / 2],
                    [DomainDimensionRatio[0] * Lm / 2, DomainDimensionRatio[1]
                        * Lm / 2, DomainDimensionRatio[2] * Lm / 2],
                ]
            )
            # Define edges of cuboid
            edges = [
                [vertices[0], vertices[1], vertices[3], vertices[2], vertices[0]],
                [vertices[4], vertices[5], vertices[7], vertices[6], vertices[4]],
                [vertices[0], vertices[4]],
                [vertices[1], vertices[5]],
                [vertices[2], vertices[6]],
                [vertices[3], vertices[7]],
            ]
            for edge in edges:
                ax.add_collection3d(Line3DCollection(
                    [edge], color="black", linewidths=1))
            # Set labels
            ax.set_box_aspect([1, 1, 1])
            ax.set_xlabel('X (m)')
            ax.set_ylabel('Y (m)')
            ax.set_zlabel('Z (m)')

            plt.show()
            return

        def MayaviVisualFlow():
            try:
                # Replace "your_cpp_executable" with the path to your C++ executable
                cpp_executable = "python3 " + current_directory + "/DFN_FLOW_VISUAL.py"
                # output = subprocess.check_output([cpp_executable], shell=True)
                # print("DFN_Gen output:")
                # print(output.decode("utf-8"))  # Decode bytes to string
                process = subprocess.Popen(cpp_executable, stdout=subprocess.PIPE,
                                           stderr=subprocess.STDOUT, universal_newlines=True, shell=True)
                # Read and print the output line by line
                for line in process.stdout:
                    print(line, end="")
                # Wait for the process to finish
                return_code = process.wait()
                # ----------
                if (return_code != 0):
                    raise subprocess.CalledProcessError(
                        return_code, cpp_executable, "Error")

            except subprocess.CalledProcessError as e:
                print("Error running DFN_Gen:", e)

                # ----------
            return

        VisualFlowWindow = Toplevel()
        VisualFlowWindow.title("Visualize flow")
        VisualFlowWindow.geometry("300x100")

        btn_visFlow = Button(
            VisualFlowWindow, text="Use Matplotlib", command=MatplotlibVisualFlow)
        btn_visFlow.pack()

        btn_visFlow2 = Button(
            VisualFlowWindow, text="Use Mayavi (require this installed)", command=MayaviVisualFlow)
        btn_visFlow2.pack()
        return
    mainWindowMesh = Toplevel()
    mainWindowMesh.title("DFN flow")
    mainWindowMesh.geometry("300x130")
    btn1 = Button(mainWindowMesh, text="Set flow parameters",
                  command=FlowParaSet)
    btn1.pack(fill=X, expand=1)

    btn2 = Button(mainWindowMesh, text="Visualize flow",
                  command=VisualizeFlowDFN)
    btn2.pack(fill=X, expand=1)
    return


def PTOptionWindow():
    def ReRunPT():
        try:
            # Replace "your_cpp_executable" with the path to your C++ executable
            cpp_executable = (
                cuDFNsys_GUI_path
                + "/DFN_PT "
                + current_directory
                + "/PTPara"
            )
            # output = subprocess.check_output(
            #    [cpp_executable], shell=True)
            # print("DFN_Gen output:")
            # print(output.decode("utf-8"))  # Decode bytes to string
            process = subprocess.Popen(
                cpp_executable, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, shell=True)
            # Read and print the output line by line
            for line in process.stdout:
                print(line, end="")
            # Wait for the process to finish
            return_code = process.wait()

            if (return_code != 0):
                raise subprocess.CalledProcessError(
                    return_code, cpp_executable, "Error")
            # ----------
            success_window = Toplevel()
            success_window.title("DFN PT")
            success_label = Label(
                success_window, text="DFN PT finished!")
            success_label.pack()
        except subprocess.CalledProcessError as e:
            print("Error running DFN_Gen:", e)
            # ----------
            success_window = Toplevel()
            success_window.title("Error")
            success_label = Label(
                success_window, text="Error happens!")
            success_label.pack()
        return

    def CalculatePeDm():
        def getDeltaT():
            Pe = entry1.get()
            LengthScale = entry2.get()
            f2 = h5py.File(current_directory + "/Class_MESH.h5")
            MeanElementArea = np.array(f2["MeanGridSize"][0])
            f2.close()
            f2 = h5py.File(current_directory + "/Class_FLOW.h5")
            MaxVelocity = np.array(f2["MaxVelocity"][0])
            MeanVelocity = np.array(f2["MeanVelocity"][0])
            f2.close()
            textaddstring = "MeanElementArea: " + \
                str(MeanElementArea) + ", MaxVelocity: " + str(MaxVelocity) + \
                ", MeanVelocity: " + \
                str(MeanVelocity) + \
                ", LengthScale you set is" + LengthScale
            textaddstring = textaddstring + "\n characterAdvectionTimeScale is: MeanElementArea^0.5 / MaxVelocity = " + \
                "{:.5e}".format(MeanElementArea ** 0.5 / MaxVelocity)
            if not Pe:
                textaddstring = textaddstring + \
                    "\n If no diffusion, then time step (delta t) is recommended to be characterAdvectionTimeScale / 10 = " + "{:.5e}".format(
                        MeanElementArea ** 0.5 / MaxVelocity/10)
            if Pe:
                Dm = float(MeanVelocity) * float(LengthScale) / float(Pe)
                textaddstring = textaddstring + \
                    "\n If you want the Peclet number (Pe) is " + Pe + \
                    ", then Dm should be: MeanVelocity * LengthScale / Pe = " + \
                    "{:.5e}".format(
                        Dm)
                textaddstring = textaddstring + "\n characterDiffusionTimeScale is: MeanElementArea / (2 * Dm) = " + str(
                    MeanElementArea / (2*Dm)) + "\n Now compare characterDiffusionTimeScale and characterAdvectionTimeScale"
                if (MeanElementArea / (2*Dm) < MeanElementArea ** 0.5 / MaxVelocity):
                    textaddstring = textaddstring + "\n characterDiffusionTimeScale is smaller, meaning the diffusion is dominant! so we better take characterDiffusionTimeScale/10 as delta t:" + \
                        "{:.5e}".format(MeanElementArea / (2*Dm)/10)
                else:
                    textaddstring = textaddstring + "\n characterAdvectionTimeScale is smaller, meaning the advection is dominant! so we better take characterAdvectionTimeScale/10 as delta t:" + \
                        "{:.5e}".format(MeanElementArea **
                                        0.5 / MaxVelocity/10)
            print("\n\n" + textaddstring + "\n\n")
            CalculatePeDmWindow.destroy()
            RunNewPT()
            return
        CalculatePeDmWindow = Toplevel()
        CalculatePeDmWindow.geometry("500x300")
        CalculatePeDmWindow.title("Get recommended time step size (delta t)")
        texts1 = r"Peclet number (Pe) you want (leave it blank for no molecular diffusion):"
        label1 = Label(CalculatePeDmWindow, text=texts1, wraplength=450)
        label1.pack()
        entry1 = Entry(CalculatePeDmWindow)
        entry1.pack()
        texts2 = r"Set the magnitude of the length scale in Pe:"
        label2 = Label(CalculatePeDmWindow, text=texts2, wraplength=450)
        label2.pack()
        entry2 = Entry(CalculatePeDmWindow)
        entry2.pack()
        button_runmore = Button(
            CalculatePeDmWindow, text="Get recommended delta t", command=getDeltaT)
        button_runmore.pack()
        return

    def RunNewPT():

        def StartRun():
            generate_or_empty_file(current_directory
                                   + "/PTPara.csv")
            if os.path.exists(current_directory + "/ParticlePositionResult"):
                shutil.rmtree(current_directory + "/ParticlePositionResult")
            if os.path.exists(current_directory + "/Dispersion_MeanSquareDisplacement.h5"):
                os.remove(current_directory +
                          "/Dispersion_MeanSquareDisplacement.h5")
            if os.path.exists(current_directory + "/EdgesSharedEle.h5"):
                os.remove(current_directory + "/EdgesSharedEle.h5")
            if os.path.exists(current_directory + "/FracturesForParticle.h5"):
                os.remove(current_directory + "/FracturesForParticle.h5")

            NumParticles = entry_widgets[0].get()
            NumTimeSteps = entry_widgets[1].get()
            MolecularDiffusion = entry_widgets[2].get()
            DeltaT = entry_widgets[3].get()
            InjectionMethod_initialcondition = entry_widgets[4].get()
            If_OutputAllPTInformationOrFPTCurve = entry_widgets[5].get()
            SpacingOfControlPlanes = entry_widgets[6].get()
            IfOutputVarianceOfDisplacementsEachStep = entry_widgets[7].get()
            IfInjectAtCustomedPlane = entry_widgets[8].get()
            CustomedPlaneInjection = entry_widgets[9].get()
            IfUseFluxWeightedOrEqualProbableMixingIntersection = entry_widgets[10].get(
            )

            if not MolecularDiffusion:
                MolecularDiffusion = "0"
            if not InjectionMethod_initialcondition:
                InjectionMethod_initialcondition = "Flux-weighted"
            if not If_OutputAllPTInformationOrFPTCurve:
                If_OutputAllPTInformationOrFPTCurve = "0"
            if not SpacingOfControlPlanes:
                SpacingOfControlPlanes = "1e7"
            if not IfOutputVarianceOfDisplacementsEachStep:
                IfOutputVarianceOfDisplacementsEachStep = "0"
            if not IfInjectAtCustomedPlane:
                IfInjectAtCustomedPlane = "0"
            if not CustomedPlaneInjection:
                CustomedPlaneInjection = "0"
            if not IfUseFluxWeightedOrEqualProbableMixingIntersection:
                IfUseFluxWeightedOrEqualProbableMixingIntersection = "1"

            NumParticles = "NumParticles," + NumParticles
            NumTimeSteps = "NumTimeSteps," + NumTimeSteps
            MolecularDiffusion = "MolecularDiffusion," + MolecularDiffusion
            DeltaT = "DeltaT," + DeltaT
            InjectionMethod_initialcondition = "InjectionMethod_initialcondition," + \
                InjectionMethod_initialcondition
            If_OutputAllPTInformationOrFPTCurve = "If_OutputAllPTInformationOrFPTCurve," + \
                If_OutputAllPTInformationOrFPTCurve
            SpacingOfControlPlanes = "SpacingOfControlPlanes," + SpacingOfControlPlanes
            IfOutputVarianceOfDisplacementsEachStep = "IfOutputVarianceOfDisplacementsEachStep," + \
                IfOutputVarianceOfDisplacementsEachStep
            IfInjectAtCustomedPlane = "IfInjectAtCustomedPlane," + IfInjectAtCustomedPlane
            CustomedPlaneInjection = "CustomedPlaneInjection," + CustomedPlaneInjection
            IfUseFluxWeightedOrEqualProbableMixingIntersection = "IfUseFluxWeightedOrEqualProbableMixingIntersection," + \
                IfUseFluxWeightedOrEqualProbableMixingIntersection

            add_line_to_file(current_directory +
                             "/PTPara.csv", NumParticles + ",")
            add_line_to_file(current_directory +
                             "/PTPara.csv", NumTimeSteps + ",")
            add_line_to_file(current_directory + "/PTPara.csv",
                             MolecularDiffusion + ",")
            add_line_to_file(current_directory + "/PTPara.csv", DeltaT + ",")
            add_line_to_file(current_directory + "/PTPara.csv",
                             InjectionMethod_initialcondition + ",")
            add_line_to_file(current_directory + "/PTPara.csv",
                             If_OutputAllPTInformationOrFPTCurve + ",")
            add_line_to_file(current_directory + "/PTPara.csv",
                             SpacingOfControlPlanes + ",")
            add_line_to_file(current_directory + "/PTPara.csv",
                             IfOutputVarianceOfDisplacementsEachStep + ",")
            add_line_to_file(current_directory + "/PTPara.csv",
                             IfInjectAtCustomedPlane + ",")
            add_line_to_file(current_directory + "/PTPara.csv",
                             CustomedPlaneInjection + ",")
            add_line_to_file(current_directory + "/PTPara.csv",
                             IfUseFluxWeightedOrEqualProbableMixingIntersection + ",")

            ReRunPT()
            ptParaWindow.destroy()
            return

        # ptParaWindow = Toplevel()

        ptParaWindow = Toplevel()
        ptParaWindow.geometry("500x800")
        ptParaWindow.title("Set PT parameters")

        texts = [
            "Number of expected particles",
            "Number of time steps (delta t)",
            "Modelcular diffusion coefficient (Dm)(leave it blank for 0):",
            "Time step size (delta t):",
            "Injection method (three options only: Flux-weighted, Resident, Point) (leave it blank for Flux-weighted):",
            "Do you want to output particles' coordinates at each step (1 for yes, 0 for first passage time only) (leave it blank for 0):",
            "Spacing of control planes (leave it blank for no control planes):",
            "Do you want to output the variance of particle displacements at each step (1 for yes) (leave it blank for 0):",
            "Do you want to specify the injection plane (1 for yes) (leave it blank for 0):",
            "The location of the injection plane (if the above entry is filled with 1, then this entry should be filled with a number which is the coordinate of a injection plane):",
            "At fracture intersections, mixing is 0 (for equiprobable mixing) or 1 (outgoing-flow-flux-weighted) (leave it blank for 1):"
        ]

        entry_widgets = []

        # Define function to create label and entry for each text
        def create_label_entry(row, text):
            label = Label(ptParaWindow, text=text, wraplength=250)
            label.grid(row=row, column=0, sticky='w', padx=5, pady=5)
            entry = Entry(ptParaWindow)
            entry.grid(row=row, column=1, padx=5, pady=5)
            entry_widgets.append(entry)  # Store the Entry widget in the list
            return entry

        # Create labels and entries in two columns
        for i, text in enumerate(texts):
            entry_sPT = create_label_entry(i, text)

        button_PT_run = Button(
            ptParaWindow, text="Run a new PT", command=StartRun)
        button_PT_run.grid(row=len(texts), columnspan=2, padx=5, pady=5)

        # label_below_button = Label(
        #    ptParaWindow, text=textaddstring, wraplength=400)
        # label_below_button.grid(
        #    row=len(texts) + 1, columnspan=2, padx=5, pady=5)
        return

    def RunMoreSteps():
        def RunMore():
            NumSteps = entry1.get()
            with open(current_directory + "/PTPara.csv", 'r') as file:
                lines = file.readlines()
            lines[1] = "NumTimeSteps, " + NumSteps + ',\n'
            with open(current_directory + "/PTPara.csv", 'w') as file:
                file.writelines(lines)
            ReRunPT()
            RunMoreStepsWindow.destroy()
            return

        RunMoreStepsWindow = Toplevel()
        RunMoreStepsWindow.geometry("300x200")
        RunMoreStepsWindow.title("Run more steps")

        text1 = r"Number of steps:"
        label1 = Label(RunMoreStepsWindow, text=text1, wraplength=250)
        label1.pack()
        entry1 = Entry(RunMoreStepsWindow)
        entry1.pack()

        button_runmore = Button(
            RunMoreStepsWindow, text="Continue", command=RunMore)
        button_runmore.pack()
        return

    def Tran2DTo3D():
        try:
            # Replace "your_cpp_executable" with the path to your C++ executable
            cpp_executable = (
                cuDFNsys_GUI_path
                + "/Transform2DH5ParticleDataTo3D  0 "
                + current_directory
                + "/DFN_MESH_VISUAL.h5"
            )
            # output = subprocess.check_output(
            #    [cpp_executable], shell=True)
            # print("DFN_Gen output:")
            # print(output.decode("utf-8"))  # Decode bytes to string
            process = subprocess.Popen(
                cpp_executable, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, shell=True)
            # Read and print the output line by line
            for line in process.stdout:
                print(line, end="")
            # Wait for the process to finish
            return_code = process.wait()

            if (return_code != 0):
                raise subprocess.CalledProcessError(
                    return_code, cpp_executable, "Error")
            # ----------
            success_window = Toplevel()
            success_window.title("DFN PT 2D to 3D")
            success_label = Label(
                success_window, text="DFN PT 2D to 3D finished!")
            success_label.pack()
        except subprocess.CalledProcessError as e:
            print("Error running Transform2DH5ParticleDataTo3D:", e)
            # ----------
            success_window = Toplevel()
            success_window.title("Error")
            success_label = Label(
                success_window, text="Error happens!")
            success_label.pack()
        return

    def VisualizePTLastStep():
        def UseMatplotlibPT():
            f = h5py.File(current_directory + "/Class_DFN.h5")
            NumFracs = np.array(f["NumFractures"][0])
            Lm = np.array(f["L"][0])
            DomainDimensionRatio = np.array(f["DomainDimensionRatio"][:])
            # print(NumFracs)
            Verts = np.zeros((NumFracs, 4, 3))
            # Faces = []

            # Create a 3D plot
            fig = plt.figure()
            ax = fig.add_subplot(111, projection="3d")
            MaxR = 0

            for i in range(0, NumFracs):
                vsd = np.array(f["Fracture_" + str(i + 1) + "/Verts3D"][:])
                vsd = np.transpose(vsd)
                Verts[i, :, :] = vsd
                R_rr = np.linalg.norm(vsd[1, :] - vsd[3, :])
                if R_rr > MaxR:
                    MaxR = R_rr
                # Faces.append([0, 1, 2, 3])
            for verts in Verts:
                poly = Poly3DCollection(
                    [verts], facecolors="blue", linewidths=1, edgecolors=(0, 0, 0, 0.1), alpha=0.1
                )
                ax.add_collection3d(poly)

            # Define coordinates of cuboid vertices
            vertices = np.array(
                [
                    [-DomainDimensionRatio[0] * Lm / 2, -DomainDimensionRatio[1]
                        * Lm / 2, -DomainDimensionRatio[2] * Lm / 2],
                    [-DomainDimensionRatio[0] * Lm / 2, -DomainDimensionRatio[1]
                        * Lm / 2, DomainDimensionRatio[2] * Lm / 2],
                    [-DomainDimensionRatio[0] * Lm / 2, DomainDimensionRatio[1]
                        * Lm / 2, -DomainDimensionRatio[2] * Lm / 2],
                    [-DomainDimensionRatio[0] * Lm / 2, DomainDimensionRatio[1]
                        * Lm / 2, DomainDimensionRatio[2] * Lm / 2],
                    [DomainDimensionRatio[0] * Lm / 2, -DomainDimensionRatio[1]
                        * Lm / 2, -DomainDimensionRatio[2] * Lm / 2],
                    [DomainDimensionRatio[0] * Lm / 2, -DomainDimensionRatio[1]
                        * Lm / 2, DomainDimensionRatio[2] * Lm / 2],
                    [DomainDimensionRatio[0] * Lm / 2, DomainDimensionRatio[1]
                        * Lm / 2, -DomainDimensionRatio[2] * Lm / 2],
                    [DomainDimensionRatio[0] * Lm / 2, DomainDimensionRatio[1]
                        * Lm / 2, DomainDimensionRatio[2] * Lm / 2],
                ]
            )
            # Define edges of cuboid
            edges = [
                [vertices[0], vertices[1], vertices[3], vertices[2], vertices[0]],
                [vertices[4], vertices[5], vertices[7], vertices[6], vertices[4]],
                [vertices[0], vertices[4]],
                [vertices[1], vertices[5]],
                [vertices[2], vertices[6]],
                [vertices[3], vertices[7]],
            ]
            for edge in edges:
                ax.add_collection3d(Line3DCollection(
                    [edge], color="black", linewidths=1))

            ax.set_xlim(
                [
                    -DomainDimensionRatio[0] * Lm / 2.0 - R_rr / 2.0,
                    DomainDimensionRatio[0] * Lm / 2.0 + R_rr / 2.0,
                ]
            )
            ax.set_ylim(
                [
                    -DomainDimensionRatio[1] * Lm / 2.0 - R_rr / 2.0,
                    DomainDimensionRatio[1] * Lm / 2.0 + R_rr / 2.0,
                ]
            )
            ax.set_zlim(
                [
                    -DomainDimensionRatio[2] * Lm / 2.0 - R_rr / 2.0,
                    DomainDimensionRatio[2] * Lm / 2.0 + R_rr / 2.0,
                ]
            )
            f.close()
            # ------------------------------PT
            f_2 = h5py.File(current_directory +
                            "/ParticlePositionResult/DispersionInfo.h5")
            N_steps = int(np.array(f_2['NumOfSteps'][0]))
            # N_particles = int(np.array(f_2['NumParticles'][0]))
            BlockNOPresent = int(np.array(f_2['BlockNOPresent'][0]))
            # SizeOfDataBlock = int(np.array(f_2['SizeOfDataBlock'][0]))
            H5name = "/ParticlePositionResult/ParticlePositionBlock" + \
                str(math.ceil(BlockNOPresent)).zfill(10) + "_3D.h5"
            f_3 = h5py.File(current_directory + H5name)
            particlePos = np.array(f_3["Step_" + str(N_steps).zfill(10)][:])
            particlePos = np.transpose(particlePos)
            ax.scatter(particlePos[:, 0], particlePos[:, 1],
                       particlePos[:, 2], c='r', marker='o', alpha=1)

            f_3.close()
            f_2.close()
            # -------------------------------
            ax.set_box_aspect([1, 1, 1])
            ax.set_xlabel("X (m)")
            ax.set_ylabel("Y (m)")
            ax.set_zlabel("Z (m)")
            ax.set_title(
                "The finall position of all particles; fractures are transparent")
            plt.show()
            return

        def UseMataviPT():
            try:
                # Replace "your_cpp_executable" with the path to your C++ executable
                cpp_executable = ("python3 " +
                                  current_directory
                                  + "/DFN_DISPERSION.py "
                                  )
                # output = subprocess.check_output(
                #    [cpp_executable], shell=True)
                # print("DFN_Gen output:")
                # print(output.decode("utf-8"))  # Decode bytes to string
                process = subprocess.Popen(
                    cpp_executable, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, shell=True)
                # Read and print the output line by line
                for line in process.stdout:
                    print(line, end="")
                # Wait for the process to finish
                return_code = process.wait()

                if (return_code != 0):
                    raise subprocess.CalledProcessError(
                        return_code, cpp_executable, "Error")
                # ----------

            except subprocess.CalledProcessError as e:
                print("Error running PT mayavi animation visualization:", e)

            return

        VisualizePTLastStepWindow = Toplevel()
        VisualizePTLastStepWindow.title("Visualize PT")
        VisualizePTLastStepWindow.geometry("300x300")

        text2 = "Only works when you set `Do you want to output particles' coordinates at each step` to 1! And the particle positions have been transformed to 3D"
        label2 = Label(VisualizePTLastStepWindow, text=text2, wraplength=250)
        label2.pack()

        btn1 = Button(VisualizePTLastStepWindow,
                      text="Use Matplotlib (just the finnal position)", command=UseMatplotlibPT)
        btn1.pack(fill=X, expand=1)

        btn1 = Button(VisualizePTLastStepWindow,
                      text="Use Mayavi (can see animation if the 2D particle displacement is tranformed to 3D)", command=UseMataviPT, wraplength=250)
        btn1.pack(fill=X, expand=1)

        text1 = "The Matlab is recommended to use to visualize PT annimation, just run `DFN_DISPERSION.m` in Matlab"
        label1 = Label(VisualizePTLastStepWindow, text=text1, wraplength=250)
        label1.pack()

        return
    mainWindowPT = Toplevel()
    mainWindowPT.title("Particle tracking")
    mainWindowPT.geometry("300x200")

    btn1 = Button(mainWindowPT, text="Run a new PT", command=CalculatePeDm)
    btn1.pack(fill=X, expand=1)

    btn2 = Button(
        mainWindowPT, text="Run more steps for an existing PT", command=RunMoreSteps)
    btn2.pack(fill=X, expand=1)

    btn4 = Button(
        mainWindowPT, text="Transform particle displacements to 3D", command=Tran2DTo3D)
    btn4.pack(fill=X, expand=1)

    btn3 = Button(mainWindowPT, text="Visualize PT",
                  command=VisualizePTLastStep)
    btn3.pack(fill=X, expand=1)

    return


def Help_window():
    stringHelp = "Please go to https://github.com/qq1012510777/cuDFNsys/blob/main/Manual/Manual.md (require VPN at mainland, China) or contact yintingchang@foxmail.com"
    print(stringHelp)
    helpwindow = Toplevel()
    helpwindow.title("Help")
    helpwindow.geometry("300x170")
    label3 = Label(helpwindow, text=stringHelp, wraplength=270)
    label3.pack()

    return


def DefineCUDFNSYS_root():
    global cuDFNsys_root_dir
    global cuDFNsys_GUI_path
    global btn1
    global btn2
    global btn3
    global btn4
    global btn4
    stringui = entry1ds.get()
    if not stringui:
        home_directory = os.path.expanduser('~')
        # print("Home directory:", home_directory)
        stringui = home_directory
        cuDFNsys_root_dir = stringui + "/cuDFNsys"
        print("cuDFNsys root directory:", cuDFNsys_root_dir)
        cuDFNsys_GUI_path = cuDFNsys_root_dir + "/GUI"
    for widget in root.winfo_children():
        widget.destroy()

    root.geometry("300x300")
    root.title("cuDFNsys")
    btn1 = Button(root, text="DFN generation", command=DFN_input_para)

    btn2 = Button(root, text="Mesh generation", command=MeshOption)

    btn3 = Button(root, text="Flow solver", command=FlowOptionWindow)

    btn4 = Button(root, text="Particle tracking", command=PTOptionWindow)

    btn5 = Button(root, text="Help", command=Help_window)
    btn1.pack(fill=X, expand=1)
    btn2.pack(fill=X, expand=1)
    btn3.pack(fill=X, expand=1)
    btn4.pack(fill=X, expand=1)
    btn5.pack(fill=X, expand=1)

    # Create a secondary window
    secondary_window = Toplevel(root)
    secondary_window.title("Output information")

    # Create a text widget for displaying terminal output
    terminal_output = Text(secondary_window, wrap=WORD)
    terminal_output.pack(fill=BOTH, expand=1)

    # Redirect stdout to the text widget
    redirect_output(terminal_output)

    root.update_idletasks()
    root_pos = [root.winfo_rootx(), root.winfo_rooty()]
    root_width, root_height = root.winfo_width(), root.winfo_height()
    secondary_window.geometry(f"+{root_pos[0] + root_width}+{root_pos[1]}")


def redirect_output(text_widget):
    class StdoutRedirector:
        def __init__(self, text_widget):
            self.text_space = text_widget

        def write(self, str):
            self.text_space.insert(END, str)
            self.text_space.see(END)  # Autoscroll to the bottom
            self.text_space.update_idletasks()  # Update the text widget

    sys.stdout = StdoutRedirector(text_widget)


# 1
root = Tk()
# 2
root.geometry("300x140")
root.title("cuDFNsys")

label3 = Label(root, text="Please specify the root directory of `cuDFNsys`(leave it blank if `cuDFNsys` is under you home directory):", wraplength=250)
label3.pack()

entry1ds = Entry(root)
entry1ds.pack()


btns3 = Button(root, text="Continue", command=DefineCUDFNSYS_root)
btns3.pack(fill=X, expand=1)


root.mainloop()
