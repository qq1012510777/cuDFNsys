from tkinter import *
import subprocess
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

def generate_or_empty_file(file_name):
    # Open the file in write mode ('w'), which truncates the file if it exists or creates a new empty file if it doesn't
    with open(file_name, "w"):
        pass  # Using pass since we don't need to perform any write operation


def add_line_to_file(file_name, line):
    with open(file_name, "a") as file:  # Open the file in append mode ('a')
        file.write(line + "\n")


current_directory = os.path.dirname(os.path.abspath(__file__))


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
                    current_directory + "/DFN_Gen " + csv_name + " " + seed_random
                )
                output = subprocess.check_output([cpp_executable], shell=True)
                print("DFN_Gen output:")
                print(output.decode("utf-8"))  # Decode bytes to string
                new_window.destroy()
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
                    cpp_executable = current_directory + "/DFN_Gen " + file_name
                    output = subprocess.check_output(
                        [cpp_executable], shell=True)
                    print("DFN_Gen output:")
                    print(output.decode("utf-8"))  # Decode bytes to string
                    new_window.destroy()
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

                    button1 = Button(
                        StochasticWindow_1, text="Generate DFN", command=(GetRandomseed)
                    )
                    button1.pack()

                    try:
                        # Replace "your_cpp_executable" with the path to your C++ executable
                        cpp_executable = (
                            current_directory
                            + "/DFN_Gen "
                            + current_directory
                            + "/StochasticFracs "
                            + randomseed_sd
                        )
                        output = subprocess.check_output(
                            [cpp_executable], shell=True)
                        print("DFN_Gen output:")
                        print(output.decode("utf-8"))  # Decode bytes to string
                        StochasticWindow_1.destroy()
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
        StochasticWindow_1.geometry("300x300")

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
                    [-DomainDimensionRatio[0] * Lm / 2, -DomainDimensionRatio[1] * Lm / 2, -DomainDimensionRatio[2] * Lm / 2],
                    [-DomainDimensionRatio[0] * Lm / 2, -DomainDimensionRatio[1] * Lm / 2, DomainDimensionRatio[2] * Lm / 2],
                    [-DomainDimensionRatio[0] * Lm / 2, DomainDimensionRatio[1] * Lm / 2, -DomainDimensionRatio[2] * Lm / 2],
                    [-DomainDimensionRatio[0] * Lm / 2, DomainDimensionRatio[1] * Lm / 2, DomainDimensionRatio[2] * Lm / 2],
                    [DomainDimensionRatio[0] * Lm / 2, -DomainDimensionRatio[1] * Lm / 2, -DomainDimensionRatio[2] * Lm / 2],
                    [DomainDimensionRatio[0] * Lm / 2, -DomainDimensionRatio[1] * Lm / 2, DomainDimensionRatio[2] * Lm / 2],
                    [DomainDimensionRatio[0] * Lm / 2, DomainDimensionRatio[1] * Lm / 2, -DomainDimensionRatio[2] * Lm / 2],
                    [DomainDimensionRatio[0] * Lm / 2, DomainDimensionRatio[1] * Lm / 2, DomainDimensionRatio[2] * Lm / 2],
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
                ax.add_collection3d(Line3DCollection([edge], color="black", linewidths=1))

            
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

            ax.set_xlabel("X (m)")
            ax.set_ylabel("Y (m)")
            ax.set_zlabel("Z (m)")

            plt.show()
            visualizeDFN_window.destroy()

        def UseMayavi():
            try:
                # Replace "your_cpp_executable" with the path to your C++ executable
                cpp_executable = "python3 " + current_directory + "/DFN_VISUAL.py"
                output = subprocess.check_output([cpp_executable], shell=True)
                # print("DFN_Gen output:")
                print(output.decode("utf-8"))  # Decode bytes to string
                visualizeDFN_window.destroy()
                # ----------

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
    new_window.title("DFN input parameters")
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


# 1
root = Tk()
# 2
root.geometry("300x300")
root.title("cuDFNsys")

btn1 = Button(root, text="DFN generation", command=DFN_input_para)
btn1.pack(fill=X, expand=1)

btn2 = Button(root, text="Mesh generation")
btn2.pack(fill=X, expand=1)

btn3 = Button(root, text="Flow solver")
btn3.pack(fill=X, expand=1)

btn4 = Button(root, text="Particle tracking")
btn4.pack(fill=X, expand=1)

btn4 = Button(root, text="Help")
btn4.pack(fill=X, expand=1)

root.mainloop()
