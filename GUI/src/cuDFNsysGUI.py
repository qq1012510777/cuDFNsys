from tkinter import *
from tkinter import filedialog
from tkinter.scrolledtext import ScrolledText
import multiprocessing

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
from pathlib import Path

import time

cuDFNsys_root_dir = ""
cuDFNsys_GUI_path = ""
cuDFNsys_root_dir = os.path.expanduser("~") + "/cuDFNsys"
cuDFNsys_GUI_path = cuDFNsys_root_dir + "/GUI"


def generate_or_empty_file(file_name):
    # Open the file in write mode ('w'), which truncates the file if it exists or creates a new empty file if it doesn't
    with open(file_name, "w"):
        pass  # Using pass since we don't need to perform any write operation


def add_line_to_file(file_name, line):
    with open(file_name, "a") as file:  # Open the file in append mode ('a')
        file.write(line + "\n")


def PrintRedString(string_to_print):
    output_text.insert(END, string_to_print + "\n", "red")
    output_text.see(END)  # Autoscroll to the bottom


def PrintGreenString(string_to_print):
    output_text.insert(END, string_to_print + "\n", "green")
    output_text.see(END)  # Autoscroll to the bottom


def RedirectOutput(text_widget):
    class StdRedirector:
        def __init__(self, text_widget, stream, redirect_output=True):
            self.text_widget = text_widget
            self.stream = stream
            self.redirect_output = redirect_output

        def write(self, message):
            if self.redirect_output:
                self.text_widget.insert(END, message)
                self.text_widget.see(END)  # Autoscroll to the bottom
                self.text_widget.update_idletasks()

        def flush(self):
            pass  # No flushing needed for GUI

    sys.stdout = StdRedirector(text_widget, sys.stdout)
    sys.stderr = StdRedirector(text_widget, sys.stderr)


def ClearAllWidget():
    for widget in root.winfo_children():
        if widget not in (filemenu, output_frame, output_text, menu):
            widget.destroy()


def OpenProject():
    currentdir = filedialog.askdirectory()
    os.chdir(currentdir)
    print("Working directory changed to:", os.getcwd())
    return


def SetcuDFNsysRoot():
    global cuDFNsys_root_dir
    global cuDFNsys_GUI_path

    def SetRoot():
        global cuDFNsys_root_dir
        global cuDFNsys_GUI_path
        cuDFNsys_root_dir = ""
        cuDFNsys_GUI_path = ""
        stringPath = entry1.get()
        if len(stringPath) != 0:
            cuDFNsys_root_dir = stringPath
        else:
            cuDFNsys_root_dir = os.path.expanduser("~")
        cuDFNsys_root_dir = cuDFNsys_root_dir + "/cuDFNsys"
        print("cuDFNsys root directory:", cuDFNsys_root_dir)
        cuDFNsys_GUI_path = cuDFNsys_root_dir + "/GUI"
        ClearAllWidget()
        return

    label1 = Label(root, width=20, height=1, text="Input cuDFNsys root:")
    label1.place(relx=0.0, rely=0.0, anchor="nw")
    entry1 = Entry(root, width=12)
    entry1.place(relx=0.2, rely=0.0, anchor="nw")

    Helpbutton = Button(
        root,
        text="Help",
        fg="Red",
        command=lambda: PrintRedString(
            "The path of where is `cuDFNsys`, e.g., /home/tingchangyin. Leave it blank for home directory"
        ),
    )
    Helpbutton.place(relx=0.35, rely=0.0, anchor="nw")

    button = Button(root, text="Set cuDFNsys root", command=(SetRoot))
    button.place(relx=0.0, rely=0.05, anchor="nw")

    return


def NewProject():
    def CreateProject():
        projectname = entry1.get()
        os.mkdir(currentdir + "/" + projectname)
        os.chdir(currentdir + "/" + projectname)
        print("Working directory changed to:", os.getcwd())
        ClearAllWidget()
        return

    # Ask the user to select a directory
    currentdir = filedialog.askdirectory()
    # If the user cancels, return without changing the working directory
    if not currentdir:
        return
    # Change the working directory
    label1 = Label(root, width=12, height=1, text="Project name:")
    label1.place(relx=0.0, rely=0.0, anchor="nw")
    entry1 = Entry(root, width=12)
    entry1.place(relx=0.15, rely=0.0, anchor="nw")

    button = Button(root, text="Create project", command=(CreateProject))
    button.place(relx=0.0, rely=0.05, anchor="nw")


def GenFractures(CSVName, randomseed):
    randomseed22 = (
        randomseed if len(
            randomseed) != 0 else "use current time as random seed"
    )
    print("Running DFN_Gen " + CSVName + ", random seed = " + str(randomseed22))
    try:
        # Replace "your_cpp_executable" with the path to your C++ executable
        cpp_executable = (
            cuDFNsys_GUI_path + "/DFN_Gen " + CSVName + " " + str(randomseed)
        )
        process = subprocess.Popen(
            cpp_executable,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            shell=True,
        )
        # Read and print the output line by line
        for line in process.stdout:
            print(line, end="")
        # Wait for the process to finish
        return_code = process.wait()

        if return_code != 0:
            raise subprocess.CalledProcessError(
                return_code, cpp_executable, "Error")
        # ----------
        PrintGreenString("DFN generated")
    except subprocess.CalledProcessError as e:
        PrintRedString("Error running DFN_Gen: " + e)
        # ----------
    return


def DeterministicDFN(If_gen=True):
    ClearAllWidget()
    global entrys1_2
    global entrys1_3
    global entrys2_1
    global entrys2_2
    global entrys2_3
    global entrys4
    global entrys5
    global entrys6_1
    global entrys6_2
    global entrys6_3
    global entrys1_1
    global entrys3

    def AddFractures(FracNumPresent, NumFrac_total):
        global entrys1_2
        global entrys1_3
        global entrys2_1
        global entrys2_2
        global entrys2_3
        global entrys4
        global entrys5
        global entrys6_1
        global entrys6_2
        global entrys6_3
        global entrys1_1
        global entrys3
        # print("FracNumPresent: ", FracNumPresent)
        if FracNumPresent > 1:

            string_frac = (
                "Fracture_"
                + str(FracNumPresent - 1)
                + ", "
                + entrys1_1.get()
                + ", "
                + entrys1_2.get()
                + ", "
                + entrys1_3.get()
                + ", "
                + entrys2_1.get()
                + ", "
                + entrys2_2.get()
                + ", "
                + entrys2_3.get()
                + ", "
            )
            string_frac = (
                string_frac
                + entrys3.get()
                + ", "
                + entrys4.get()
                + ", "
                + entrys5.get()
                + ","
            )
            string_frac = (
                string_frac
                + entrys6_1.get()
                + ", "
                + entrys6_2.get()
                + ", "
                + entrys6_3.get()
                + ","
            )
            add_line_to_file("./DeterministicDFN.csv", string_frac)
            ClearAllWidget()
        if FracNumPresent > NumFrac_total:
            if If_gen:
                GenFractures("./DeterministicDFN", "")
            return

        if FracNumPresent == 1:
            Lm = float(entry1.get())

            DomainRatio_x = 1 if len(
                entry2_1.get()) == 0 else float(entry2_1.get())
            DomainRatio_y = 1 if len(
                entry2_2.get()) == 0 else float(entry2_2.get())
            DomainRatio_y = (
                DomainRatio_y if DomainRatio_x == 1 else DomainRatio_y / DomainRatio_x
            )
            DomainRatio_z = 1 if len(
                entry2_3.get()) == 0 else float(entry2_3.get())
            DomainRatio_z = (
                DomainRatio_z if DomainRatio_x == 1 else DomainRatio_z / DomainRatio_x
            )

            PercolationDirection = (
                2 if (len(entry3.get()) == 0) else float(entry3.get())
            )
            NumFractures = 1 if (len(entry4.get()) == 0) else int(entry4.get())

            generate_or_empty_file("./DeterministicDFN.csv")
            add_line_to_file("./DeterministicDFN.csv", "IfStochastic,0,")

            add_line_to_file("./DeterministicDFN.csv",
                             "DomainSizeX, " + str(Lm) + ",")
            add_line_to_file(
                "./DeterministicDFN.csv",
                "DomainDimensionRatio,"
                + str(DomainRatio_x)
                + ","
                + str(DomainRatio_y)
                + ","
                + str(DomainRatio_z)
                + ",",
            )
            add_line_to_file(
                "./DeterministicDFN.csv",
                "Percolation_direction," + str(PercolationDirection) + ",",
            )
            add_line_to_file(
                "./DeterministicDFN.csv", "NumFractures," +
                str(NumFractures) + ","
            )
            NumFrac_total = NumFractures
            ClearAllWidget()

        # Fracture normal vector
        relateX = 0
        relateY = 0.0
        labels1 = Label(root, width=20, height=1, text="Normal vector:")
        labels1.place(relx=relateX, rely=relateY, anchor="nw")
        entrys1_1 = Entry(root, width=2)
        entrys1_1.place(relx=relateX + 0.2, rely=relateY, anchor="nw")
        entrys1_2 = Entry(root, width=2)
        entrys1_2.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
        entrys1_3 = Entry(root, width=2)
        entrys1_3.place(relx=relateX + 0.3, rely=relateY, anchor="nw")
        Helpbuttons1 = Button(
            root,
            text="Help",
            fg="Red",
            height=1,
            command=lambda: PrintRedString(
                "Input three values to make the normal vector of the fracture"
            ),
        )
        Helpbuttons1.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

        # Fracture normal vector
        relateX = 0
        relateY = 0.05
        labels2 = Label(root, width=20, height=1, text="Center:")
        labels2.place(relx=relateX, rely=relateY, anchor="nw")
        entrys2_1 = Entry(root, width=2)
        entrys2_1.place(relx=relateX + 0.2, rely=relateY, anchor="nw")
        entrys2_2 = Entry(root, width=2)
        entrys2_2.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
        entrys2_3 = Entry(root, width=2)
        entrys2_3.place(relx=relateX + 0.3, rely=relateY, anchor="nw")
        Helpbuttons2 = Button(
            root,
            text="Help",
            fg="Red",
            height=1,
            command=lambda: PrintRedString(
                "Input three values to make the center of the fracture"
            ),
        )
        Helpbuttons2.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

        # size
        relateX = 0
        relateY = 0.1
        labels3 = Label(root, width=20, height=1, text="Size:")
        labels3.place(relx=relateX, rely=relateY, anchor="nw")
        entrys3 = Entry(root, width=12)
        entrys3.place(relx=relateX + 0.2, rely=relateY, anchor="nw")
        Helpbuttons3 = Button(
            root,
            text="Help",
            fg="Red",
            height=1,
            command=lambda: PrintRedString(
                "The radius of the circumscribed circle of the square fracture"
            ),
        )
        Helpbuttons3.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

        # beta
        relateX = 0
        relateY = 0.15
        labels4 = Label(root, width=20, height=1, text="Beta:")
        labels4.place(relx=relateX, rely=relateY, anchor="nw")
        entrys4 = Entry(root, width=12)
        entrys4.place(relx=relateX + 0.2, rely=relateY, anchor="nw")
        Helpbuttons4 = Button(
            root,
            text="Help",
            fg="Red",
            height=1,
            command=lambda: PrintRedString(
                "both Beta and Gamma relate to the apeter and conductivity of fractures. Note that the conductivity: k_f = b_f^3 / 12 =  [(Gamma * R ^ Beta)] ^ 3 / 12)"
            ),
        )
        Helpbuttons4.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

        # gamma
        relateX = 0
        relateY = 0.2
        labels5 = Label(root, width=20, height=1, text="Gamma:")
        labels5.place(relx=relateX, rely=relateY, anchor="nw")
        entrys5 = Entry(root, width=12)
        entrys5.place(relx=relateX + 0.2, rely=relateY, anchor="nw")
        Helpbuttons5 = Button(
            root,
            text="Help",
            fg="Red",
            height=1,
            command=lambda: PrintRedString(
                "both Beta and Gamma relate to the apeter and conductivity of fractures. Note that the conductivity: k_f = b_f^3 / 12 =  [(Gamma * R ^ Beta)] ^ 3 / 12)"
            ),
        )
        Helpbuttons5.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

        # specify one vertex
        relateX = 0
        relateY = 0.25
        labels6 = Label(root, width=20, height=1, text="One vertex:")
        labels6.place(relx=relateX, rely=relateY, anchor="nw")
        entrys6_1 = Entry(root, width=2)
        entrys6_1.place(relx=relateX + 0.2, rely=relateY, anchor="nw")
        entrys6_2 = Entry(root, width=2)
        entrys6_2.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
        entrys6_3 = Entry(root, width=2)
        entrys6_3.place(relx=relateX + 0.3, rely=relateY, anchor="nw")
        Helpbuttons6 = Button(
            root,
            text="Help",
            fg="Red",
            height=1,
            command=lambda: PrintRedString(
                "Specify one vertex of the square fracture. Leave them blank for random selection of this vertex."
            ),
        )
        Helpbuttons6.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

        # button
        button_text_add_frac = (
            ("Add fracture " + str(FracNumPresent))
            if FracNumPresent < NumFrac_total
            else ("Generate")
        )

        buttonAdd_fractures = Button(
            root,
            text=button_text_add_frac,
            command=(lambda: AddFractures(FracNumPresent, NumFrac_total)),
        )
        buttonAdd_fractures.place(relx=relateX, rely=0.3, anchor="nw")
        FracNumPresent = FracNumPresent + 1

    # domain size x
    label1 = Label(root, width=20, height=1, text="DomainsizeX:")
    label1.place(relx=0.0, rely=0.0, anchor="nw")
    entry1 = Entry(root, width=12)
    entry1.place(relx=0.2, rely=0.0, anchor="nw")
    Helpbutton1 = Button(
        root,
        text="Help",
        fg="Red",
        height=1,
        command=lambda: PrintRedString("Domain size in the x direction"),
    )
    Helpbutton1.place(relx=0.35, rely=0.0, anchor="nw")

    # domain dimension ratio
    relateX = 0
    relateY = 0.05
    label2 = Label(root, width=20, height=1, text="Domain dimension ratio:")
    label2.place(relx=relateX, rely=relateY, anchor="nw")
    entry2_1 = Entry(root, width=2)
    entry2_1.place(relx=relateX + 0.2, rely=relateY, anchor="nw")
    entry2_2 = Entry(root, width=2)
    entry2_2.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
    entry2_3 = Entry(root, width=2)
    entry2_3.place(relx=relateX + 0.3, rely=relateY, anchor="nw")
    Helpbutton2 = Button(
        root,
        text="Help",
        fg="Red",
        height=1,
        command=lambda: PrintRedString(
            "Leave them blank for 1, 1, 1, i.e., a cubic domain. If you input 1, 1, 2 here and the DomainsizeX is 30, then the domain size is, 30 m x 30 m x 60 m"
        ),
    )
    Helpbutton2.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

    # percolation direction
    relateX = 0
    relateY = 0.1
    label3 = Label(root, width=20, height=1, text="Percolation direction:")
    label3.place(relx=relateX, rely=relateY, anchor="nw")
    entry3 = Entry(root, width=3)
    entry3.place(relx=relateX + 0.2, rely=relateY, anchor="nw")
    Helpbutton3 = Button(
        root,
        text="Help",
        fg="Red",
        height=1,
        command=lambda: PrintRedString(
            "Which direction is the percolative direction. 0 for x, 1 for y, 2 for z. Leave it blank for 2."
        ),
    )
    Helpbutton3.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

    # Number of fractures
    relateX = 0
    relateY = 0.15
    label4 = Label(root, width=20, height=1, text="Number of fractures:")
    label4.place(relx=relateX, rely=relateY, anchor="nw")
    entry4 = Entry(root, width=3)
    entry4.place(relx=relateX + 0.2, rely=relateY, anchor="nw")
    Helpbutton4 = Button(
        root,
        text="Help",
        fg="Red",
        height=1,
        command=lambda: PrintRedString("Leave it blank to be 1"),
    )
    Helpbutton4.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

    # Add fractures
    FracNumPresent = 1
    global NumFrac_total
    NumFrac_total = 1
    buttonAdd_fractures = Button(
        root,
        text="Add fractures",
        command=(lambda: AddFractures(FracNumPresent, NumFrac_total)),
    )
    buttonAdd_fractures.place(relx=0, rely=0.2, anchor="nw")


def StochasticDFN(If_gen=True):
    ClearAllWidget()
    global entry_t_1
    global entry_t_2
    global entrys_t_3_1
    global entrys_t_3_2
    global entrys_t_3_3
    global entry_t_4
    global entrys_t_5_1
    global entrys_t_5_2
    global entrys_t_5_3
    global entrys_t_5_4
    global entry_t_6
    global entry_t_7
    RandomSeed_stochastic_genDFN = ""

    def AddFractureGroups(GroupPresent, NumberGroups):
        global NumFractureEachGroup_str
        global KappaValues_str
        global MeanOrientationOfFisherDistribution_str
        global ModeOfSizeDistribution_str
        global SizeDistributionParameters_str
        global Beta_str
        global Gamma_str
        global entry_t_1
        global entry_t_2
        global entrys_t_3_1
        global entrys_t_3_2
        global entrys_t_3_3
        global entry_t_4
        global entrys_t_5_1
        global entrys_t_5_2
        global entrys_t_5_3
        global entrys_t_5_4
        global entry_t_6
        global entry_t_7
        global RandomSeed_stochastic_genDFN

        if GroupPresent > 1:
            NumFractureEachGroup_str = NumFractureEachGroup_str + entry_t_1.get() + ", "
            KappaValues_str = KappaValues_str + entry_t_2.get() + ", "
            MeanOrientationOfFisherDistribution_str = (
                MeanOrientationOfFisherDistribution_str
                + entrys_t_3_1.get()
                + ", "
                + entrys_t_3_2.get()
                + ", "
                + entrys_t_3_3.get()
                + ", "
            )
            ModeOfSizeDistribution_str = (
                ModeOfSizeDistribution_str + entry_t_4.get() + ", "
            )
            SizeDistributionParameters_str = (
                SizeDistributionParameters_str
                + entrys_t_5_1.get()
                + ", "
                + entrys_t_5_2.get()
                + ", "
                + entrys_t_5_3.get()
                + ", "
                + entrys_t_5_4.get()
                + ", "
            )
            Beta_str = Beta_str + entry_t_6.get() + ", "
            Gamma_str = Gamma_str + entry_t_7.get() + ", "
            ClearAllWidget()
        if GroupPresent > NumberGroups:
            add_line_to_file("./StochasticDFN.csv", NumFractureEachGroup_str)
            add_line_to_file("./StochasticDFN.csv", KappaValues_str)
            add_line_to_file(
                "./StochasticDFN.csv", MeanOrientationOfFisherDistribution_str
            )
            add_line_to_file("./StochasticDFN.csv", ModeOfSizeDistribution_str)
            add_line_to_file("./StochasticDFN.csv",
                             SizeDistributionParameters_str)
            add_line_to_file("./StochasticDFN.csv", Beta_str)
            add_line_to_file("./StochasticDFN.csv", Gamma_str)
            if If_gen:
                GenFractures("./StochasticDFN", RandomSeed_stochastic_genDFN)
            return

        if GroupPresent == 1:
            Lm = float(entry1.get())

            DomainRatio_x = 1 if len(
                entry2_1.get()) == 0 else float(entry2_1.get())
            DomainRatio_y = 1 if len(
                entry2_2.get()) == 0 else float(entry2_2.get())
            DomainRatio_y = (
                DomainRatio_y if DomainRatio_x == 1 else DomainRatio_y / DomainRatio_x
            )
            DomainRatio_z = 1 if len(
                entry2_3.get()) == 0 else float(entry2_3.get())
            DomainRatio_z = (
                DomainRatio_z if DomainRatio_x == 1 else DomainRatio_z / DomainRatio_x
            )

            PercolationDirection = (
                2 if (len(entry3.get()) == 0) else float(entry3.get())
            )
            NumGroups_lo = 1 if (len(entry4.get()) == 0) else int(entry4.get())

            RandomSeed_stochastic_genDFN = entry5.get()

            generate_or_empty_file("./StochasticDFN.csv")
            add_line_to_file("./StochasticDFN.csv", "IfStochastic,1,")

            add_line_to_file("./StochasticDFN.csv",
                             "DomainSizeX, " + str(Lm) + ",")
            add_line_to_file(
                "./StochasticDFN.csv",
                "DomainDimensionRatio,"
                + str(DomainRatio_x)
                + ","
                + str(DomainRatio_y)
                + ","
                + str(DomainRatio_z)
                + ",",
            )
            add_line_to_file(
                "./StochasticDFN.csv",
                "Percolation_direction," + str(PercolationDirection) + ",",
            )
            add_line_to_file(
                "./StochasticDFN.csv", "NumFractures," +
                str(NumGroups_lo) + ","
            )
            NumberGroups = NumGroups_lo
            ClearAllWidget()

        # Num fractures each group
        relateX = 0
        relateY = 0.0
        label_t_1 = Label(root, width=20, height=1,
                          text="Number of fractures:")
        label_t_1.place(relx=relateX, rely=relateY, anchor="nw")
        entry_t_1 = Entry(root, width=2)
        entry_t_1.place(relx=relateX + 0.2, rely=relateY, anchor="nw")
        Helpbutton_t_1 = Button(
            root,
            text="Help",
            fg="Red",
            height=1,
            command=lambda: PrintRedString(
                "Input the number of fractures in this group"
            ),
        )
        Helpbutton_t_1.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

        # kappa (fisher)
        relateX = 0
        relateY = 0.05
        label_t_2 = Label(root, width=20, height=1, text="Kappa:")
        label_t_2.place(relx=relateX, rely=relateY, anchor="nw")
        entry_t_2 = Entry(root, width=2)
        entry_t_2.place(relx=relateX + 0.2, rely=relateY, anchor="nw")
        Helpbutton_t_2 = Button(
            root,
            text="Help",
            fg="Red",
            height=1,
            command=lambda: PrintRedString(
                "Input the Fisher constant of the Fisher distribution for orientation dispersivity"
            ),
        )
        Helpbutton_t_2.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

        # Mean orientation
        relateX = 0
        relateY = 0.1
        labels_t_3 = Label(root, width=20, height=1, text="Mean orientation:")
        labels_t_3.place(relx=relateX, rely=relateY, anchor="nw")
        entrys_t_3_1 = Entry(root, width=2)
        entrys_t_3_1.place(relx=relateX + 0.2, rely=relateY, anchor="nw")
        entrys_t_3_2 = Entry(root, width=2)
        entrys_t_3_2.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
        entrys_t_3_3 = Entry(root, width=2)
        entrys_t_3_3.place(relx=relateX + 0.3, rely=relateY, anchor="nw")
        Helpbuttons_t_3 = Button(
            root,
            text="Help",
            fg="Red",
            height=1,
            command=lambda: PrintRedString(
                "Input three values to make the mean orientation of the group"
            ),
        )
        Helpbuttons_t_3.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

        # Mode of size distribution
        relateX = 0
        relateY = 0.15
        label_t_4 = Label(root, width=20, height=1, text="Size mode:")
        label_t_4.place(relx=relateX, rely=relateY, anchor="nw")
        entry_t_4 = Entry(root, width=2)
        entry_t_4.place(relx=relateX + 0.2, rely=relateY, anchor="nw")
        Helpbutton_t_4 = Button(
            root,
            text="Help",
            fg="Red",
            height=1,
            command=lambda: PrintRedString(
                "Mode of fracture size distribution (0: power-law, 1: lognormal, 2: uniform, 3: mono-sized)."
            ),
        )
        Helpbutton_t_4.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

        # SizeDistributionParameters
        relateX = 0
        relateY = 0.25
        labels_t_5 = Label(
            root, width=25, height=1, text="Size distribution parameter:"
        )
        labels_t_5.place(relx=relateX, rely=relateY, anchor="nw")
        entrys_t_5_1 = Entry(root, width=2)
        entrys_t_5_1.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
        entrys_t_5_2 = Entry(root, width=2)
        entrys_t_5_2.place(relx=relateX + 0.3, rely=relateY, anchor="nw")
        entrys_t_5_3 = Entry(root, width=2)
        entrys_t_5_3.place(relx=relateX + 0.35, rely=relateY, anchor="nw")
        entrys_t_5_4 = Entry(root, width=2)
        entrys_t_5_4.place(relx=relateX + 0.4, rely=relateY, anchor="nw")
        Helpbuttons_t_5 = Button(
            root,
            text="Help",
            fg="Red",
            height=1,
            command=lambda: PrintRedString(
                "Parameters of size distribution (note that there must be four values separatered by commas, e.g., if you input 2 in the above entry for uniform distribution, you can input 1, 10, 0, 0, meaning minimum and maximum are 1 and 10, more details are in Manual)"
            ),
        )
        Helpbuttons_t_5.place(relx=relateX + 0.45, rely=relateY, anchor="nw")

        # beta
        relateX = 0.45
        relateY = 0.0
        label_t_6 = Label(root, width=10, height=1, text="Beta:")
        label_t_6.place(relx=relateX, rely=relateY, anchor="nw")
        entry_t_6 = Entry(root, width=2)
        entry_t_6.place(relx=relateX + 0.2, rely=relateY, anchor="nw")
        Helpbutton_t_6 = Button(
            root,
            text="Help",
            fg="Red",
            height=1,
            command=lambda: PrintRedString(
                "both Beta and Gamma relate to the apeter and conductivity of fractures. Note that the conductivity: k_f = b_f^3 / 12 =  [(Gamma * R ^ Beta)] ^ 3 / 12)"
            ),
        )
        Helpbutton_t_6.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

        # gamma
        relateX = 0.45
        relateY = 0.05
        label_t_7 = Label(root, width=10, height=1, text="Gamma:")
        label_t_7.place(relx=relateX, rely=relateY, anchor="nw")
        entry_t_7 = Entry(root, width=2)
        entry_t_7.place(relx=relateX + 0.2, rely=relateY, anchor="nw")
        Helpbutton_t_7 = Button(
            root,
            text="Help",
            fg="Red",
            height=1,
            command=lambda: PrintRedString(
                "both Beta and Gamma relate to the apeter and conductivity of fractures. Note that the conductivity: k_f = b_f^3 / 12 =  [(Gamma * R ^ Beta)] ^ 3 / 12)"
            ),
        )
        Helpbutton_t_7.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

        # button
        relateX = 0.7
        button_text_add_group = (
            ("Add group " + str(GroupPresent))
            if GroupPresent < NumberGroups
            else ("Generate")
        )

        buttonAdd_fractures = Button(
            root,
            text=button_text_add_group,
            command=(lambda: AddFractureGroups(GroupPresent, NumberGroups)),
        )
        buttonAdd_fractures.place(relx=relateX, rely=0.25, anchor="nw")
        GroupPresent = GroupPresent + 1

    # domain size x
    label1 = Label(root, width=20, height=1, text="DomainsizeX:")
    label1.place(relx=0.0, rely=0.0, anchor="nw")
    entry1 = Entry(root, width=12)
    entry1.place(relx=0.2, rely=0.0, anchor="nw")
    Helpbutton1 = Button(
        root,
        text="Help",
        fg="Red",
        height=1,
        command=lambda: PrintRedString("Domain size in the x direction"),
    )
    Helpbutton1.place(relx=0.35, rely=0.0, anchor="nw")

    # domain dimension ratio
    relateX = 0
    relateY = 0.05
    label2 = Label(root, width=20, height=1, text="Domain dimension ratio:")
    label2.place(relx=relateX, rely=relateY, anchor="nw")
    entry2_1 = Entry(root, width=2)
    entry2_1.place(relx=relateX + 0.2, rely=relateY, anchor="nw")
    entry2_2 = Entry(root, width=2)
    entry2_2.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
    entry2_3 = Entry(root, width=2)
    entry2_3.place(relx=relateX + 0.3, rely=relateY, anchor="nw")
    Helpbutton2 = Button(
        root,
        text="Help",
        fg="Red",
        height=1,
        command=lambda: PrintRedString(
            "Leave them blank for 1, 1, 1, i.e., a cubic domain. If you input 1, 1, 2 here and the DomainsizeX is 30, then the domain size is, 30 m x 30 m x 60 m"
        ),
    )
    Helpbutton2.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

    # percolation direction
    relateX = 0
    relateY = 0.1
    label3 = Label(root, width=20, height=1, text="Percolation direction:")
    label3.place(relx=relateX, rely=relateY, anchor="nw")
    entry3 = Entry(root, width=3)
    entry3.place(relx=relateX + 0.2, rely=relateY, anchor="nw")
    Helpbutton3 = Button(
        root,
        text="Help",
        fg="Red",
        height=1,
        command=lambda: PrintRedString(
            "Which direction is the percolative direction. 0 for x, 1 for y, 2 for z. Leave it blank for 2."
        ),
    )
    Helpbutton3.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

    # Number of fracture groups
    relateX = 0
    relateY = 0.15
    label4 = Label(root, width=25, height=1, text="Number of fracture groups:")
    label4.place(relx=relateX, rely=relateY, anchor="nw")
    entry4 = Entry(root, width=3)
    entry4.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
    Helpbutton4 = Button(
        root,
        text="Help",
        fg="Red",
        height=1,
        command=lambda: PrintRedString(
            "Number of facture groups. Leave it blank to be 1"
        ),
    )
    Helpbutton4.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

    # random seed
    relateX = 0
    relateY = 0.2
    label5 = Label(root, width=25, height=1, text="Random seed:")
    label5.place(relx=relateX, rely=relateY, anchor="nw")
    entry5 = Entry(root, width=3)
    entry5.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
    Helpbutton5 = Button(
        root,
        text="Help",
        fg="Red",
        height=1,
        command=lambda: PrintRedString(
            "Leave it blank for using the current time as random seed"
        ),
    )
    Helpbutton5.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

    # Add fracture groups
    global NumberGroups
    global GroupPresent
    global NumFractureEachGroup_str
    global KappaValues_str
    global MeanOrientationOfFisherDistribution_str
    global ModeOfSizeDistribution_str
    global SizeDistributionParameters_str
    global Beta_str
    global Gamma_str

    NumFractureEachGroup_str = "NumFractureEachGroup" + ", "
    KappaValues_str = "KappaValues" + ", "
    MeanOrientationOfFisherDistribution_str = (
        "MeanOrientationOfFisherDistribution" + ", "
    )
    ModeOfSizeDistribution_str = "ModeOfSizeDistribution" + ", "
    SizeDistributionParameters_str = "SizeDistributionParameters" + ", "
    Beta_str = "Beta" + ", "
    Gamma_str = "Gamma" + ", "

    NumberGroups = 1
    GroupPresent = 1
    buttonAdd_fracture_groups = Button(
        root,
        text="Add fracture groups",
        command=(lambda: AddFractureGroups(GroupPresent, NumberGroups)),
    )
    buttonAdd_fracture_groups.place(relx=0, rely=0.25, anchor="nw")
    return RandomSeed_stochastic_genDFN


def UseExistingCSV():
    ClearAllWidget()

    def GenFractures_UseExistingCSV():
        GenFractures(entry2.get(), entry1.get())
        ClearAllWidget()
        return

    # rand dom seed
    relateX = 0
    relateY = 0.0
    label1 = Label(root, width=25, height=1, text="Random seed:")
    label1.place(relx=relateX, rely=relateY, anchor="nw")
    entry1 = Entry(root, width=3)
    entry1.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
    Helpbutton1 = Button(
        root,
        text="Help",
        fg="Red",
        height=1,
        command=lambda: PrintRedString(
            "Leave it blank for using the current time as random seed"
        ),
    )
    Helpbutton1.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

    # rand dom seed
    relateX = 0
    relateY = 0.05
    label2 = Label(root, width=25, height=1, text="CSV name:")
    label2.place(relx=relateX, rely=relateY, anchor="nw")
    entry2 = Entry(root, width=10)
    entry2.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
    Helpbutton2 = Button(
        root,
        text="Help",
        fg="Red",
        height=1,
        command=lambda: PrintRedString(
            "Input the csv name in the working directory without the suffix .csv"
        ),
    )
    Helpbutton2.place(relx=relateX + 0.45, rely=relateY, anchor="nw")

    buttonGen = Button(root, text="Generate",
                       command=GenFractures_UseExistingCSV)
    buttonGen.place(relx=0, rely=0.25, anchor="nw")

    return


def DFNvisualUseMatPlotLib():
    f = h5py.File("./Class_DFN.h5")
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
            [
                -DomainDimensionRatio[0] * Lm / 2,
                -DomainDimensionRatio[1] * Lm / 2,
                -DomainDimensionRatio[2] * Lm / 2,
            ],
            [
                -DomainDimensionRatio[0] * Lm / 2,
                -DomainDimensionRatio[1] * Lm / 2,
                DomainDimensionRatio[2] * Lm / 2,
            ],
            [
                -DomainDimensionRatio[0] * Lm / 2,
                DomainDimensionRatio[1] * Lm / 2,
                -DomainDimensionRatio[2] * Lm / 2,
            ],
            [
                -DomainDimensionRatio[0] * Lm / 2,
                DomainDimensionRatio[1] * Lm / 2,
                DomainDimensionRatio[2] * Lm / 2,
            ],
            [
                DomainDimensionRatio[0] * Lm / 2,
                -DomainDimensionRatio[1] * Lm / 2,
                -DomainDimensionRatio[2] * Lm / 2,
            ],
            [
                DomainDimensionRatio[0] * Lm / 2,
                -DomainDimensionRatio[1] * Lm / 2,
                DomainDimensionRatio[2] * Lm / 2,
            ],
            [
                DomainDimensionRatio[0] * Lm / 2,
                DomainDimensionRatio[1] * Lm / 2,
                -DomainDimensionRatio[2] * Lm / 2,
            ],
            [
                DomainDimensionRatio[0] * Lm / 2,
                DomainDimensionRatio[1] * Lm / 2,
                DomainDimensionRatio[2] * Lm / 2,
            ],
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


def DFNvisualUseMayavi():
    try:
        # Replace "your_cpp_executable" with the path to your C++ executable
        cpp_executable = "python3 ./DFN_VISUAL.py"
        # output = subprocess.check_output([cpp_executable], shell=True)
        # print("DFN_Gen output:")
        # print(output.decode("utf-8"))  # Decode bytes to string
        process = subprocess.Popen(
            cpp_executable,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            shell=True,
        )
        # Read and print the output line by line
        for line in process.stdout:
            print(line, end="")
        # Wait for the process to finish
        return_code = process.wait()
        PrintGreenString("Mayavi visualization of the DFN")
        # ----------
        if return_code != 0:
            raise subprocess.CalledProcessError(
                return_code, cpp_executable, "Error")
    except subprocess.CalledProcessError as e:
        PrintRedString("Error:")
        print("Error running DFN_Gen:", e)
        # ----------
    return


def Deterministic_and_tochastic_DFN():
    ClearAllWidget()

    def Gen_mixed():
        randomseed_x = entry1.get()
        randomseed_x = "-1" if len(randomseed_x) == 0 else randomseed_x

        string_h = "current time" if randomseed_x == "-1" else randomseed_x
        print(
            "Running DFN_Gen ./StochasticDFN ./DeterministicDFN, random seed = "
            + string_h
        )
        try:
            # Replace "your_cpp_executable" with the path to your C++ executable
            cpp_executable = (
                cuDFNsys_GUI_path
                + "/DFN_Gen ./StochasticDFN ./DeterministicDFN "
                + str(randomseed_x)
            )
            process = subprocess.Popen(
                cpp_executable,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
                shell=True,
            )
            # Read and print the output line by line
            for line in process.stdout:
                print(line, end="")
            # Wait for the process to finish
            return_code = process.wait()

            if return_code != 0:
                raise subprocess.CalledProcessError(
                    return_code, cpp_executable, "Error"
                )
            # ----------
            PrintGreenString(
                "DFN generated (both stochastic and deterministic fractures)"
            )
        except subprocess.CalledProcessError as e:
            PrintRedString("Error running DFN_Gen: " + e)
            # ----------
        ClearAllWidget()
        return

    label1 = Label(root, width=12, height=1, text="Random seed:")
    label1.place(relx=0.0, rely=0.0, anchor="nw")
    entry1 = Entry(root, width=12)
    entry1.place(relx=0.15, rely=0.0, anchor="nw")

    button = Button(root, text="Generate", command=(Gen_mixed))
    button.place(relx=0.0, rely=0.05, anchor="nw")

    return


def DFNMesh(IfGen=True):
    ClearAllWidget()

    def GenMesh():
        mixS = entry1.get()
        maxS = entry2.get()
        generate_or_empty_file("./MeshPara.csv")
        add_line_to_file("./MeshPara.csv",
                         "ExpectedMinimimGridSize, " + str(mixS) + ",")
        add_line_to_file("./MeshPara.csv",
                         "ExpectedMaximimGridSize, " + str(maxS) + ",")
        ClearAllWidget()
        if IfGen:
            try:
                cpp_executable = (
                    cuDFNsys_GUI_path + "/DFN_Mesh ./MeshPara"
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
                if (return_code != 0):
                    raise subprocess.CalledProcessError(
                        return_code, cpp_executable, "Error")
                # ----------
                PrintGreenString("DFN mesh generated")
            except subprocess.CalledProcessError as e:
                PrintRedString("Error:")
                print("Error running DFN_Gen:", e)

    # minmum grid size
    relateX = 0
    relateY = 0.0
    label1 = Label(root, width=25, height=1,
                   text="Minimum expected grid size:")
    label1.place(relx=relateX, rely=relateY, anchor="nw")
    entry1 = Entry(root, width=3)
    entry1.place(relx=relateX + 0.25, rely=relateY, anchor="nw")

    # maximum grid size
    relateX = 0
    relateY = 0.05
    label2 = Label(root, width=25, height=1,
                   text="Maximum expected grid size:")
    label2.place(relx=relateX, rely=relateY, anchor="nw")
    entry2 = Entry(root, width=3)
    entry2.place(relx=relateX + 0.25, rely=relateY, anchor="nw")

    button = Button(root, text="Generate", command=(GenMesh))
    button.place(relx=0.0, rely=0.1, anchor="nw")

    return


def VisualizeMeshMatplotlib():
    f = h5py.File("./Class_MESH.h5")
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
    f2 = h5py.File("./Class_DFN.h5")
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


def VisualizeMeshMayavi():
    try:
        # Replace "your_cpp_executable" with the path to your C++ executable
        cpp_executable = "python3 ./DFN_MESH_VISUAL.py"
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
    return


def FlowSlover(IfGen=True):
    ClearAllWidget()

    def FlowSim():
        generate_or_empty_file("./FlowPara.csv")

        add_line_to_file("./FlowPara.csv", "MuOverRhoG, " + str(entry1.get()) +
                         ",")
        add_line_to_file("./FlowPara.csv", "InletHead, " + str(entry2.get()) + ","
                         )
        add_line_to_file("./FlowPara.csv", "OutletHead, " + str(entry3.get()) +
                         ",")
        ClearAllWidget()
        if IfGen:
            try:
                cpp_executable = (
                    cuDFNsys_GUI_path + "/DFN_Flow ./FlowPara"
                )
                process = subprocess.Popen(cpp_executable, stdout=subprocess.PIPE,
                                           stderr=subprocess.STDOUT, universal_newlines=True, shell=True)
                for line in process.stdout:
                    print(line, end="")
                # Wait for the process to finish
                return_code = process.wait()
                if (return_code != 0):
                    raise subprocess.CalledProcessError(
                        return_code, cpp_executable, "Error")
                # ----------
                PrintGreenString("Flow simulation finished")
            except subprocess.CalledProcessError as e:
                PrintRedString("Error")
                print("Error running DFN_Flow:", e)
            return
    # para mu
    relateX = 0
    relateY = 0.0
    label1 = Label(root, width=20, height=1, text=r"Mu / (rho * g):")
    label1.place(relx=relateX, rely=relateY, anchor="nw")
    entry1 = Entry(root, width=3)
    entry1.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
    Helpbutton = Button(
        root,
        text="Help",
        fg="Red",
        command=lambda: PrintRedString(
            "input the result of fluid viscosity over (density times gravity)"
        ),
    )
    Helpbutton.place(relx=0.35, rely=0.0, anchor="nw")

    # maximum grid size
    relateX = 0
    relateY = 0.05
    label2 = Label(root, width=25, height=1, text="Inlet head [L]:")
    label2.place(relx=relateX, rely=relateY, anchor="nw")
    entry2 = Entry(root, width=3)
    entry2.place(relx=relateX + 0.25, rely=relateY, anchor="nw")

    # maximum grid size
    relateX = 0
    relateY = 0.1
    label3 = Label(root, width=25, height=1, text="Outlet head [L]:")
    label3.place(relx=relateX, rely=relateY, anchor="nw")
    entry3 = Entry(root, width=3)
    entry3.place(relx=relateX + 0.25, rely=relateY, anchor="nw")

    button = Button(root, text="Solve", command=(FlowSim))
    button.place(relx=0.0, rely=0.15, anchor="nw")
    return


def VisualizeFlowMatplotlib():
    f = h5py.File("./Class_MESH.h5")
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
    f3 = h5py.File("./Class_FLOW.h5")
    PressureElements = np.array(f3["PressureEles"][:])
    f3.close()
    surf = ax.plot_trisurf(Points[:, 0], Points[:, 1],
                           Points[:, 2], triangles=Triangles-1, cmap='viridis', shade=False)
    # Set colors based on PressureElements values
    surf.set_array(PressureElements.ravel())
    fig.colorbar(surf, ax=ax, label='Hydraulic head')
    f2 = h5py.File("./Class_DFN.h5")
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


def VisualizeFlowMayavi():
    try:
        # Replace "your_cpp_executable" with the path to your C++ executable
        cpp_executable = "python3 " + "./DFN_FLOW_VISUAL.py"
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
        print("Error running VisualizeFlowMayavi:", e)
        return


def RunNewPT():
    ClearAllWidget()

    def WriteCSV():
        if os.path.exists("." + "/ParticlePositionResult"):
            shutil.rmtree("." + "/ParticlePositionResult")
        if os.path.exists("." + "/Dispersion_MeanSquareDisplacement.h5"):
            os.remove("." +
                      "/Dispersion_MeanSquareDisplacement.h5")
        if os.path.exists("." + "/EdgesSharedEle.h5"):
            os.remove("." + "/EdgesSharedEle.h5")
        if os.path.exists("." + "/FracturesForParticle.h5"):
            os.remove("." + "/FracturesForParticle.h5")

        generate_or_empty_file("./PTPara.csv")
        add_line_to_file("./PTPara.csv", "NumParticles, " + entry1.get() + ",")
        add_line_to_file("./PTPara.csv", "NumTimeSteps, " + entry2.get() + ",")
        add_line_to_file(
            "./PTPara.csv", "MolecularDiffusion, " + entry2_dm.get() + ",")
        add_line_to_file("./PTPara.csv", "DeltaT, " + entry2_s.get() + ",")
        Injection_methos = ""
        if (len(entry3.get()) == 0):
            Injection_methos = "Flux-weighted"
        elif (entry3.get() == "0"):
            Injection_methos = "Flux-weighted"
        elif (entry3.get() == "1"):
            Injection_methos = "Resident"
        elif (entry3.get() == "2"):
            Injection_methos = "Point"
        add_line_to_file(
            "./PTPara.csv", "InjectionMethod_initialcondition," + Injection_methos + ",")
        if_output_all = ""
        if (len(entry4.get()) == 0):
            if_output_all = "1"
        else:
            if_output_all = entry4.get()
        add_line_to_file(
            "./PTPara.csv", "If_OutputAllPTInformationOrFPTCurve, " + if_output_all + ",")

        spacing_of_control_plane = "1e6" if (
            len(entry5.get()) == 0) else entry5.get()
        add_line_to_file(
            "./PTPara.csv", "SpacingOfControlPlanes, " + spacing_of_control_plane + ",")

        if_particle_variance = "1" if (
            len(entry6.get()) == 0) else entry6.get()
        add_line_to_file(
            "./PTPara.csv", "IfOutputVarianceOfDisplacementsEachStep, " + if_particle_variance + ",")

        if_change_injection_plane = "0" if (
            len(entry7.get()) == 0) else entry7.get()
        add_line_to_file(
            "./PTPara.csv", "IfInjectAtCustomedPlane, " + if_change_injection_plane + ",")

        customPlane = entry8.get()
        add_line_to_file(
            "./PTPara.csv", "CustomedPlaneInjection, " + customPlane + ",")

        mixingrule = "1" if (len(entry9.get()) == 0) else entry9.get()
        add_line_to_file(
            "./PTPara.csv", "IfUseFluxWeightedOrEqualProbableMixingIntersection, " + mixingrule + ",")

        RunMoreSteps()
        ClearAllWidget()

        return
    # NP
    relateX = 0
    relateY = 0.0
    label1 = Label(root, width=25, height=1,
                   text="NP:")
    label1.place(relx=relateX, rely=relateY, anchor="nw")
    entry1 = Entry(root, width=3)
    entry1.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
    Helpbutton1 = Button(
        root,
        text="Help",
        fg="Red",
        command=lambda: PrintRedString(
            "Number of particles expected to be injected:"
        ),
    )
    Helpbutton1.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

    # Dm
    relateX = 0
    relateY = 0.05
    label2_dm = Label(root, width=25, height=1,
                      text="Dm:")
    label2_dm.place(relx=relateX, rely=relateY, anchor="nw")
    entry2_dm = Entry(root, width=3)
    entry2_dm.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
    Helpbutton2_dm = Button(
        root,
        text="Help",
        fg="Red",
        command=lambda: PrintRedString(
            "Molecular diffusion coefficient"
        ),
    )
    Helpbutton2_dm.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

    # delta T
    relateX = 0
    relateY = 0.1
    label2_s = Label(root, width=25, height=1,
                     text="Delta T:")
    label2_s.place(relx=relateX, rely=relateY, anchor="nw")
    entry2_s = Entry(root, width=3)
    entry2_s.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
    Helpbutton2_s = Button(
        root,
        text="Help",
        fg="Red",
        command=lambda: PrintRedString(
            "The size of time steps:"
        ),
    )
    Helpbutton2_s.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

    # Injection method
    relateX = 0
    relateY = 0.15
    label3 = Label(root, width=25, height=1,
                   text="Injection method:")
    label3.place(relx=relateX, rely=relateY, anchor="nw")
    entry3 = Entry(root, width=3)
    entry3.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
    Helpbutton3 = Button(
        root,
        text="Help",
        fg="Red",
        command=lambda: PrintRedString(
            "Injection method (three options only: 0 for Flux-weighted, 1 for Resident, 2 for Point) (leave it blank for 0 - Flux-weighted)"
        ),
    )
    Helpbutton3.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

    # Particle coordinates
    relateX = 0
    relateY = 0.20
    label4 = Label(root, width=25, height=1,
                   text="Output particle coodinates:")
    label4.place(relx=relateX, rely=relateY, anchor="nw")
    entry4 = Entry(root, width=3)
    entry4.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
    Helpbutton4 = Button(
        root,
        text="Help",
        fg="Red",
        command=lambda: PrintRedString(
            "Do you want to output particles' coordinates at each step (1 for yes, 0 for first passage time only) (leave it blank for 1)"
        ),
    )
    Helpbutton4.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

    # spacing of control planes
    relateX = 0
    relateY = 0.25
    label5 = Label(root, width=25, height=1,
                   text="Spacing of control plane:")
    label5.place(relx=relateX, rely=relateY, anchor="nw")
    entry5 = Entry(root, width=3)
    entry5.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
    Helpbutton5 = Button(
        root,
        text="Help",
        fg="Red",
        command=lambda: PrintRedString(
            "Spacing of control planes (leave it blank for no control planes). Control plane means a plane where if particles cross this plane, then the arrival time is recorded"
        ),
    )
    Helpbutton5.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

    # variance displacements
    relateX = 0
    relateY = 0.30
    label6 = Label(root, width=25, height=1,
                   text="Particle displacement variance:")
    label6.place(relx=relateX, rely=relateY, anchor="nw")
    entry6 = Entry(root, width=3)
    entry6.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
    Helpbutton6 = Button(
        root,
        text="Help",
        fg="Red",
        command=lambda: PrintRedString(
            "Do you want to output the variance of particle displacements at each step (1 for yes) (leave it blank for 1)"
        ),
    )
    Helpbutton6.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

    # If change injection plane
    relateX = 0.45
    relateY = 0.0
    label7 = Label(root, width=25, height=1,
                   text="If change injection plane:")
    label7.place(relx=relateX, rely=relateY, anchor="nw")
    entry7 = Entry(root, width=3)
    entry7.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
    Helpbutton7 = Button(
        root,
        text="Help",
        fg="Red",
        command=lambda: PrintRedString(
            "Do you want to specify the injection plane (1 for yes) (leave it blank for 0)"
        ),
    )
    Helpbutton7.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

    # injection plane
    relateX = 0.45
    relateY = 0.05
    label8 = Label(root, width=25, height=1,
                   text="Injection plane:")
    label8.place(relx=relateX, rely=relateY, anchor="nw")
    entry8 = Entry(root, width=3)
    entry8.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
    Helpbutton8 = Button(
        root,
        text="Help",
        fg="Red",
        command=lambda: PrintRedString(
            "The location of the injection plane (if the above entry is filled with 1, then this entry should be filled with a number which is the coordinate of a injection plane)"
        ),
    )
    Helpbutton8.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

    # injection plane
    relateX = 0.45
    relateY = 0.1
    label9 = Label(root, width=25, height=1,
                   text="Mixing at intersections:")
    label9.place(relx=relateX, rely=relateY, anchor="nw")
    entry9 = Entry(root, width=3)
    entry9.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
    Helpbutton9 = Button(
        root,
        text="Help",
        fg="Red",
        command=lambda: PrintRedString(
            "At fracture intersections, mixing is 0 (for equiprobable mixing) or 1 (outgoing-flow-flux-weighted) (leave it blank for 1)"
        ),
    )
    Helpbutton9.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

    # NT
    relateX = 0.45
    relateY = 0.15
    label2 = Label(root, width=25, height=1,
                   text="NT:")
    label2.place(relx=relateX, rely=relateY, anchor="nw")
    entry2 = Entry(root, width=3)
    entry2.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
    Helpbutton2 = Button(
        root,
        text="Help",
        fg="Red",
        command=lambda: PrintRedString(
            "Number of time steps:"
        ),
    )
    Helpbutton2.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

    buttonsk = Button(root, text="PT start", command=(WriteCSV))
    buttonsk.place(relx=relateX + 0.15, rely=0.21, anchor="nw")

    return


def GetRecommendedDmAndDealtT():
    ClearAllWidget()

    def calculateDm_deltaT():

        Pe = -1 if len(entry1.get()) == 0 else float(entry1.get())
        LengthScale = float(entry2.get())
        f2 = h5py.File("./Class_MESH.h5")
        MeanElementArea = np.array(f2["MeanGridSize"][0])
        f2.close()
        f2 = h5py.File("./Class_FLOW.h5")
        MaxVelocity = np.array(f2["MaxVelocity"][0])
        MeanVelocity = np.array(f2["MeanVelocity"][0])
        f2.close()

        print("Pe: ", Pe)
        print("LengthScale: ", LengthScale)
        print("MeanElementArea: ", MeanElementArea)
        print("MaxVelocity: ", MaxVelocity)
        print("MeanVelocity: ", MeanVelocity)

        factor_s = 10
        ChateristicAdvectionTimeScale = MeanElementArea ** 0.5 / MaxVelocity

        print("ChateristicAdvectionTimeScale is MeanElementArea^0.5 / MaxVelocity:")
        PrintGreenString(str(ChateristicAdvectionTimeScale))

        if (Pe == -1):
            print("no diffusion exists, so the recommended time step size is characterAdvectionTimeScale / " + str(factor_s))
            PrintGreenString(str(ChateristicAdvectionTimeScale / factor_s))
        else:
            Dm = MeanVelocity * LengthScale / Pe
            characterDiffusionTimeScale = MeanElementArea / (2 * Dm)
            PrintGreenString("Dm: " + str(Dm))
            print("characterDiffusionTimeScale is MeanElementArea / (2 * Dm):")
            PrintGreenString(str(characterDiffusionTimeScale))
            if (characterDiffusionTimeScale < ChateristicAdvectionTimeScale):
                print("Diffusion is stronger!")
                print(
                    "so the recommended time step size is characterDiffusionTimeScale / " + str(factor_s) + ": ")
                PrintGreenString(str(characterDiffusionTimeScale/10))
            else:
                print("Advection is stronger!")
                print(
                    "so the recommended time step size is ChateristicAdvectionTimeScale / " + str(factor_s) + ": ")
                PrintGreenString(str(ChateristicAdvectionTimeScale/10))
        ClearAllWidget()
        RunNewPT()
        return
    # minmum grid size
    relateX = 0
    relateY = 0.0
    label1 = Label(root, width=25, height=1,
                   text="Pe:")
    label1.place(relx=relateX, rely=relateY, anchor="nw")
    entry1 = Entry(root, width=3)
    entry1.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
    Helpbutton1 = Button(
        root,
        text="Help",
        fg="Red",
        command=lambda: PrintRedString(
            "Peclet number (Pe) you want (leave it blank for no molecular diffusion):"
        ),
    )
    Helpbutton1.place(relx=0.35, rely=0.0, anchor="nw")

    # maximum grid size
    relateX = 0
    relateY = 0.05
    label2 = Label(root, width=25, height=1,
                   text="Length scale:")
    label2.place(relx=relateX, rely=relateY, anchor="nw")
    entry2 = Entry(root, width=3)
    entry2.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
    Helpbutton2 = Button(
        root,
        text="Help",
        fg="Red",
        command=lambda: PrintRedString(
            "Set the magnitude of the length scale in Pe:"
        ),
    )
    Helpbutton2.place(relx=relateX+0.35, rely=relateY, anchor="nw")

    button = Button(root, text="Generate", command=(calculateDm_deltaT))
    button.place(relx=0.0, rely=0.1, anchor="nw")

    return


def RunMoreSteps():
    try:
        # Replace "your_cpp_executable" with the path to your C++ executable
        cpp_executable = (
            cuDFNsys_GUI_path
            + "/DFN_PT "
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
        PrintGreenString("PT finished")
    except subprocess.CalledProcessError as e:
        PrintRedString("Error:")
        print("Error running DFN_Gen:", e)
        # ----------
    return


def ReRun():
    def ReRunII():
        NumSteps = entry1.get()
        with open("./PTPara.csv", 'r') as file:
            lines = file.readlines()
        lines[1] = "NumTimeSteps, " + NumSteps + ',\n'
        with open("./PTPara.csv", 'w') as file:
            file.writelines(lines)
        RunMoreSteps()
        return
    # NP
    relateX = 0
    relateY = 0.0
    label1 = Label(root, width=25, height=1,
                   text="NT:")
    label1.place(relx=relateX, rely=relateY, anchor="nw")
    entry1 = Entry(root, width=3)
    entry1.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
    Helpbutton1 = Button(
        root,
        text="Help",
        fg="Red",
        command=lambda: PrintRedString(
            "How many steps you want to continue to run?"
        ),
    )
    Helpbutton1.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

    button = Button(root, text="Solve", command=(ReRunII))
    button.place(relx=relateX, rely=0.05, anchor="nw")

    return


def Get3DDisplacement():
    try:
        # Replace "your_cpp_executable" with the path to your C++ executable
        cpp_executable = (
            cuDFNsys_GUI_path
            + "/Transform2DH5ParticleDataTo3D  0 "
            + "./DFN_MESH_VISUAL.h5"
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
        PrintGreenString("Get3DDisplacement finishied")
    except subprocess.CalledProcessError as e:
        PrintRedString("Error:")
        print("Error running Transform2DH5ParticleDataTo3D:", e)
    return


def VisualizePTMatplotlib():
    f = h5py.File("./Class_DFN.h5")
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
    f_2 = h5py.File("./ParticlePositionResult/DispersionInfo.h5")
    N_steps = int(np.array(f_2['NumOfSteps'][0]))
    # N_particles = int(np.array(f_2['NumParticles'][0]))
    BlockNOPresent = int(np.array(f_2['BlockNOPresent'][0]))
    # SizeOfDataBlock = int(np.array(f_2['SizeOfDataBlock'][0]))
    H5name = "./ParticlePositionResult/ParticlePositionBlock" + \
        str(math.ceil(BlockNOPresent)).zfill(10) + "_3D.h5"
    f_3 = h5py.File(H5name)
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


def VisualizePTMayavi():
    try:
        # Replace "your_cpp_executable" with the path to your C++ executable
        cpp_executable = ("python3 "
                          + "./DFN_DISPERSION.py "
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


def SetDFNParameters():
    def CallStochasticDFN():
        StochasticDFN(False)
        return

    relateX = 0.25
    relateY = 0.05
    label1 = Label(root, height=1,
                   text="MC iterations for DFN and flow simulation!", fg="red", font=("Helvetica", 14))
    label1.place(relx=relateX, rely=relateY, anchor="nw")

    button = Button(root, text="Continue", command=(CallStochasticDFN))
    button.place(relx=0.5, rely=0.1, anchor="nw")

    return


def RunMCIterations():

    ClearAllWidget()

    def StartRunMC():
        NUMMC = 100 if (len(entry3.get()) == 0) else int(entry3.get())
        NumProcessor = 5 if (len(entry4.get()) == 0) else int(entry4.get())
        DevGPUInx = 0 if (len(entry5.get()) == 0) else int(entry5.get())
        PrintGreenString("Using GPU card No. " + str(DevGPUInx))
        ClearAllWidget()

        if (NumProcessor < 1):
            PrintRedString("Number of CPU threads must be >= 1")
            return
        if (NumProcessor == 1):
            PrintGreenString("Number of iterations: " +
                             str(NUMMC) + ", Serial MC iteration ...")
        if (NumProcessor > 1):
            PrintGreenString("Number of iterations: " + str(NUMMC) +
                             ", NumProcessor = " + str(NumProcessor) + ", Parallel MC iteration ... ")

        def countdown(seconds, string):
            while seconds > 0:
                print(".")
                time.sleep(1)
                seconds -= 1
                if (seconds == 1):
                    PrintGreenString(string)
        countdown(3, "Start MC iterations!")

        try:
            # Replace "your_cpp_executable" with the path to your C++ executable
            cpp_executable = (
                cuDFNsys_GUI_path
                + "/MCIterations "
                + cuDFNsys_GUI_path + " " +
                str(NUMMC) + " " + str(NumProcessor) + " " + str(DevGPUInx)
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
            PrintGreenString("MC iterations finished")
        except subprocess.CalledProcessError as e:
            PrintRedString("Error:")
            print("Error running MCIterations:", e)
        # ----------
        PrintGreenString("MC finished!")
        return
    # m
    relateX = 0
    relateY = 0.0
    label3 = Label(root, width=25, height=1, text="Number of MC iterations:")
    label3.place(relx=relateX, rely=relateY, anchor="nw")
    entry3 = Entry(root, width=3)
    entry3.place(relx=relateX + 0.25, rely=relateY, anchor="nw")

    #
    relateX = 0
    relateY = 0.05
    label4 = Label(root, width=25, height=1,
                   text="Number of CPU threads (NCT):")
    label4.place(relx=relateX, rely=relateY, anchor="nw")
    entry4 = Entry(root, width=3)
    entry4.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
    Helpbutton4 = Button(
        root,
        text="Help",
        fg="Red",
        command=lambda: PrintRedString(
            "There will be NCT tasks running at the same time, depending on one GPU card"
        ),
    )
    Helpbutton4.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

    #
    relateX = 0
    relateY = 0.1
    label5 = Label(root, width=25, height=1,
                   text="The index of GPU card:")
    label5.place(relx=relateX, rely=relateY, anchor="nw")
    entry5 = Entry(root, width=3)
    entry5.place(relx=relateX + 0.25, rely=relateY, anchor="nw")
    Helpbutton5 = Button(
        root,
        text="Help",
        fg="Red",
        command=lambda: PrintRedString(
            "Leave it blank for the first GPU card, i.e., No, 0. Or input 1, 2, or 3 if you have multiple GPU cards!"
        ),
    )
    Helpbutton5.place(relx=relateX + 0.35, rely=relateY, anchor="nw")

    button = Button(root, text="Start", command=StartRunMC)
    button.place(relx=relateX + 0.1, rely=relateY + 0.15, anchor="nw")
    return


def CleanAllLogs():
    ClearAllWidget()
    for root, dirs, files in os.walk("./"):
        for file in files:
            if "log_MC_" in file:
                file_path = os.path.join(root, file)
                os.remove(file_path)
                print(f"Deleted file: {file_path}")
    return


def CleanAllMCData():
    ClearAllWidget()

    def CleanStart():
        ClearAllWidget()
        for root, dirs, files in os.walk("./", topdown=False):
            for dir_name in dirs:
                if "MC_" in dir_name:
                    dir_path = os.path.join(root, dir_name)
                    shutil.rmtree(dir_path)
                    print(f"Deleted directory: {dir_path}")
        return

    def exitCleanMCData():
        ClearAllWidget()
        return
    relateX = 0
    relateY = 0.0
    label4 = Label(root, width=70, height=1,
                   text="Are you sure to delete all MC directories in the current directory?", fg="red")
    label4.place(relx=relateX, rely=relateY, anchor="nw")

    button = Button(root, text="Yes", command=CleanStart)
    button.place(relx=relateX, rely=relateY + 0.05, anchor="nw")

    button2 = Button(root, text="No", command=exitCleanMCData)
    button2.place(relx=relateX, rely=relateY + 0.1, anchor="nw")

    return


def Contact():
    PrintGreenString(
        "Please go to https://github.com/qq1012510777/cuDFNsys/blob/main/Manual/Manual.md (require VPN at mainland, China) or contact yintingchang@foxmail.com")
    return


# --------------------root window
root = Tk()
root.geometry("800x700")
root.title("cuDFNsys v1.0.1")
menu = Menu(root)
root.config(menu=menu)
filemenu = Menu(menu)

# Create a frame to contain the output text
output_frame = Frame(root, height=200)
# Create a Text widget for displaying the output
output_text = ScrolledText(output_frame, wrap=WORD)
output_text.pack(fill=BOTH, expand=True)
output_text.tag_config("red", foreground="red")
output_text.tag_config("green", foreground="green")

# place frame
output_frame.pack(side=BOTTOM, fill=X, expand=False)

# Redirect stdout to the Text widget
RedirectOutput(output_text)

print("Default root directory of cuDFNsys is " + cuDFNsys_root_dir)

menu.add_cascade(label="File", menu=filemenu)
filemenu.add_command(label="New project", command=NewProject)
filemenu.add_command(label="Open project", command=OpenProject)
filemenu.add_command(label="Set cuDFNsys root", command=SetcuDFNsysRoot)

# ---------------DFN gen
DFNgenMenu = Menu(menu)
menu.add_cascade(label="DFN generation", menu=DFNgenMenu)
DFNgenMenu.add_command(label="Deterministic DFN", command=DeterministicDFN)
DFNgenMenu.add_command(label="Stochastic DFN", command=StochasticDFN)
DFNgenMenu.add_command(label="Use existing .csv", command=UseExistingCSV)

MixedDFN_submenu = Menu(DFNgenMenu)
MixedDFN_submenu.add_command(
    label="Set stochastic fractures", command=(lambda: StochasticDFN(False))
)
MixedDFN_submenu.add_command(
    label="Set deterministic fractures", command=(lambda: DeterministicDFN(False))
)
MixedDFN_submenu.add_command(
    label="Generate", command=Deterministic_and_tochastic_DFN)
DFNgenMenu.add_cascade(
    label="Deterministic and Stochastic DFN", menu=MixedDFN_submenu)

# Create a submenu for the VisualizeDFNMenu
visualize_submenu = Menu(DFNgenMenu)
# Add commands or sub-cascades to the submenu
visualize_submenu.add_command(
    label="Use Matplotlib", command=DFNvisualUseMatPlotLib)
visualize_submenu.add_command(
    label="Use Mayavi2 (if installed)", command=DFNvisualUseMayavi
)
# Add the submenu as a cascade to the VisualizeDFNMenu
DFNgenMenu.add_cascade(label="Visualize DFN", menu=visualize_submenu)

# ---------------Mesh
DFNmeshMenu = Menu(menu)
menu.add_cascade(label="Mesh", menu=DFNmeshMenu)
DFNmeshMenu.add_command(label="Mesh", command=DFNMesh)

MeshVisual_submenu = Menu(DFNmeshMenu)
MeshVisual_submenu.add_command(
    label="Use Matplotlib", command=VisualizeMeshMatplotlib)
MeshVisual_submenu.add_command(
    label="Use Mayavi2 (if installed)", command=VisualizeMeshMayavi)
DFNmeshMenu.add_cascade(label="Visualize mesh", menu=MeshVisual_submenu)

# ---------------Flow
DFNFlowMenu = Menu(menu)
menu.add_cascade(label="Flow", menu=DFNFlowMenu)
DFNFlowMenu.add_command(label="Set flow parameter", command=FlowSlover)

FlowVisual_submenu = Menu(DFNFlowMenu)
FlowVisual_submenu.add_command(
    label="Use Matplotlib", command=VisualizeFlowMatplotlib)
FlowVisual_submenu.add_command(
    label="Use Mayavi2 (if installed)", command=VisualizeFlowMayavi)
DFNFlowMenu.add_cascade(label="Visualize mesh", menu=FlowVisual_submenu)

# ---------------PT
DFNPTMenu = Menu(menu)
menu.add_cascade(label="Particle tracking", menu=DFNPTMenu)
DFNPTMenu.add_command(label="Run new PT", command=GetRecommendedDmAndDealtT)
DFNPTMenu.add_command(
    label="Run more steps for existing PT", command=ReRun)
DFNPTMenu.add_command(label="Get 3D particle displacements",
                      command=Get3DDisplacement)

PTVisual_submenu = Menu(DFNPTMenu)
PTVisual_submenu.add_command(
    label="Use Matplotlib", command=VisualizePTMatplotlib)
PTVisual_submenu.add_command(
    label="Use Mayavi2 (if installed)", command=VisualizePTMayavi)
DFNPTMenu.add_cascade(label="Visualize mesh", menu=PTVisual_submenu)

# ---------------MonteCarlo simulation
DFNMC_SIM_Menu = Menu(menu)
menu.add_cascade(label="Monte Carlo iterations", menu=DFNMC_SIM_Menu)
DFNMC_SIM_Menu.add_command(label="Set DFN inputs", command=SetDFNParameters)
DFNMC_SIM_Menu.add_command(label="Set Mesh inputs",
                           command=(lambda: DFNMesh(False)))
DFNMC_SIM_Menu.add_command(label="Set Flow inputs",
                           command=(lambda: FlowSlover(False)))
DFNMC_SIM_Menu.add_command(label="Run MC iterations", command=RunMCIterations)

DFNMC_SIM_Menu.add_command(label="Clean all log files", command=CleanAllLogs)

DFNMC_SIM_Menu.add_command(
    label="Clean all data files", command=CleanAllMCData)

# ---------------Help
HelpMenu = Menu(menu)
menu.add_cascade(label="Help", menu=HelpMenu)
HelpMenu.add_command(label="Contact", command=Contact)

mainloop()
