#!/usr/local/bin/python3

import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib.patches as mpatches

import matplotlib.cm as cm

#------------------------------------------------------------------------------------#

def file_len(fname):
  with open(fname) as f:
    for i, l in enumerate(f):
      pass
    return i + 1

#------------------------------------------------------------------------------------#
def plot_phase_space(plot_file, npt):

    n_lines = file_len(plot_file)
    print("reading ", plot_file, ", found n =", n_lines," lines")

    plt.ylabel(r'$(a_tQ)^2$')
    plt.xlabel(r'$a_tEcm$')

    f = open(plot_file, 'r')

    data_in = f.read()

    data = data_in.split("\n\n")

    for idx,dat in enumerate(data, 1):
        if(npt == 1 and idx == 2): plot_dat = dat
        if(npt == 3 and idx == 3): plot_dat = dat  
        else:
            continue
            

#------------------------------------------------------------------------------------#

    plot_data = plot_dat.split("\n")

    for elem in plot_data:

        this_line = (elem.rstrip('\n')).split()

        if ( len(this_line) == 0):
            continue
        if (this_line[0]=="#"):
            #commented out lines begin with "# "
            continue

        ecm = float(this_line[0])
        xerr = float(this_line[1])
        qsq = float(this_line[2])
        yerr  = float(this_line[3])

	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

        plt.errorbar( ecm, qsq, yerr=yerr,xerr=xerr, markersize=10,fmt='o',color='black', mfc='w',mec='black',elinewidth=0.5, capsize=6, mew=1,zorder =10)


    plt.title(r"Phase Space")

    # leg = mpatches.Patch(color='black', label='qsq')

    # plt.legend(handles=[leg])


    fontP = FontProperties()
    fontP.set_size('small')
    plt.tight_layout()

#------------------------------------------------------------------------------------#

# ----- begin main program ----- #

if (len(sys.argv)!=3):
    print("usage: ", str(sys.argv[0]), "<plot_file> <npt>")
    exit(1)

plot_file = str(sys.argv[1])
npt       = int(sys.argv[2])

fig = plt.figure(figsize=(20, 10))
plot_phase_space(plot_file, npt)

fig.patch.set_facecolor('white')
plt.savefig("phase_space.pdf", bbox_inches='tight')
plt.show()
exit(0)

# End of plotter
