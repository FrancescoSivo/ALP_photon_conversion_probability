#Code for representing the behaviour of the probability with respect to one parameter
import math
import matplotlib.pyplot as plt
#import matplotlib.colors as colors

def readMatrix(fn):
    """Function to input the data from the file

    Args:
        fn (str): The name of the file to be read

    Returns:
        matrix: The matrix of the data
    """
    matrix = []
    with open(fn) as f:                 # a with block will auto close your file after the statements within it
        for line in f:
            line = line.strip()         # strip off any trailing whitespace(including '\n')
            matrix.append(line.split()) 
    return matrix

def progress_bar(progress, total):
    """Function to print a progress bar to the console

    Args:
        progress (float): The current progress of the task
        total (float): The total number of steps in the task
    """
    percent=100*(progress / float(total))
    bar = '█'*int(percent) + '-'*(100-int(percent))   # "█" = alt+219
    print(f"\r|{bar}| {percent:.2f}%", end="\r")
    if percent >= 99.99:
        print(f"\r|{bar}| {percent:.2f}%", end="\r")

#Inporting the data from the file and elaboring them
tab = readMatrix("1D_plot.txt")
M = 0.0
m = 1.0
MB = 0.0
ty = float(tab[0][0])
tb = 0
if len(tab[0])>1:
    tb = int(tab[0][1])
if ty==3:
    scale = float(tab[0][1])
xvec = []
p = []
b = []
N = len(tab)
n = len(tab[1])

progress_bar(0, N-1)
for i in range(1,N):
    xvec.append(float(tab[i][0]))
    p.append(float(tab[i][1]))
    if n>2:
        b.append(float(tab[i][2]))
        if float(tab[i][2]) >= MB:
            MB = float(tab[i][2])
    if float(tab[i][1]) >= M:
        M = float(tab[i][1])
    if float(tab[i][1]) <= m:
        m = float(tab[i][1])
    progress_bar(i, N-1)
print("\n")

order = math.floor(math.log(M, 10))
if n>2:
    orderB = math.floor(math.log(MB, 10))
else:
    orderB = 0

for i in range(N-1):
    p[i]*=10**(-order)
if orderB!=0 and n>2:
    for i in range(N-1):
        b[i]*=10**(-orderB)

font = {'family': 'serif',
        'color':  'Black',
        'weight': 'normal',
        'size': 14,
        }
font_blue = {'family': 'serif',
        'color':  'Blue',
        'weight': 'normal',
        'size': 14,
        }
font_red = {'family': 'serif',
        'color':  'Red',
        'weight': 'normal',
        'size': 14,
        }

pstring = "P" + " (x10$^" + str(-order) + "$)"

if n>2:
    if orderB!=0:
        if tb==1:
            bstring = "B" + " (x10$^" + str(-orderB) + "$ $\mu$G)"
        elif tb==2:
            bstring = "B$_T$" + " (x10$^" + str(-orderB) + "$ $\mu$G)"   
    else:
        if tb==1:
            bstring = "B ($\mu$G)"
        elif tb==2:
            bstring = "B$_T$ ($\mu$G)" 
if n==2:
    plt.plot(xvec, p, 'b-')
    if ty==1:
        plt.xlabel("z (kpc)", fontdict=font)
    elif ty==2:
        plt.xlabel("ma (neV)", fontdict=font)
    elif ty==3:
        if scale==0:
            plt.xlabel("E (eV)", fontdict=font)
        elif scale==1:
            plt.xlabel("E (keV)", fontdict=font)
        elif scale==2:
            plt.xlabel("E (MeV)", fontdict=font)
        elif scale==3:
            plt.xlabel("E (GeV)", fontdict=font)
        elif scale==4:
            plt.xlabel("E (TeV)", fontdict=font)
        elif scale==5:
            plt.xlabel("E (PeV)", fontdict=font)
    elif ty==4:
        plt.xlabel("g$_{a\gamma}$ (GeV$^{-1}$)", fontdict=font)
    if ty!=1:
        plt.xscale('log')
    if(m!=0):
        if(math.log(M,10)-math.log(m,10))>2.0:
            plt.yscale('log')
    plt.ylabel(pstring, fontdict=font)
    plt.grid()
    plt.savefig("1D_plot.pdf")
else:
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(xvec, p, 'b-')
    ax2.plot(xvec, b, 'r-')
    ax1.set_xlabel('z (kpc)', fontdict=font)
    ax1.set_ylabel(pstring, fontdict=font_blue)
    ax2.set_ylabel(bstring, fontdict=font_red)
    ax1.grid()
    plt.savefig("1D_plot.pdf")

print("Process completed.")
#end of code