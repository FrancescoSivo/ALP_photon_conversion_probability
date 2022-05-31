#Code with the purpose of producing a Mollweide projection of the data
#importing the necessary libraries
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors


def progress_bar(progress, total):
    """Function to print a progress bar to the console

    Args:
        progress (float): The current progress of the task
        total (float): The total number of steps in the task
    """
    percent=100*(progress / float(total))
    bar = '█'*int(percent) + '-'*(100-int(percent))   # "█" = alt+219
    print( f"\r|{bar}| {percent:.2f}%", end="\r")
    if percent >= 99.99:
        print( f"\r|{bar}| {percent:.2f}%", end="\r")

def readMatrix(fn):
    """Function to input the data from the file

    Args:
        fn (str): The name of the file to be read
    """
    matrix = []
    with open(fn) as f:                 # a with block will auto close your file after the statements within it
        for line in f:
            line = line.strip()         # strip off any trailing whitespace(including '\n')
            matrix.append(line.split()) 
    return matrix

#Inporting the data from the file and elaboring them
tab = readMatrix("skymap.txt")
skymap = []
scalesetting = int(tab[0][0])
minprob = float(tab[1][0])
minprobscale = float(tab[1][1])
minprob*=10**-minprobscale
maxprob = float(tab[2][0])
maxprobscale = float(tab[2][1])
maxprob*=10**-maxprobscale
for i in range(3,len(tab)):
    skymap.append(tab[i])
print("skymap opened correctly")
print("Generating the Mollweide Plot...")
divmappa = int(math.sqrt(len(skymap)))-1
for row in skymap:
    row[0]=float(row[0])
    row[1]=float(row[1])
    row[2]=float(row[2])

# Rewriting the data to the correct format for the skymap
arr1 = np.zeros([divmappa + 1,divmappa + 1])
lon = np.zeros(divmappa + 1)
lat = np.zeros(divmappa + 1)
i = 0
while i<=divmappa:
    j = 0
    while j<=divmappa:
        arr1[j][divmappa-i] = skymap[i*(divmappa + 1) + j][0]
        if j==0:
            lon[i]=skymap[i*(divmappa + 1) + j][1]*math.pi/180
        if i==0:
            lat[j]=skymap[i*(divmappa + 1) + j][2]*math.pi/180
        progress_bar((i*(divmappa + 1) + j) , (divmappa+1)**2 )
        j += 1
    i += 1
progress_bar(1 , 1)
print('\n')

# Generating the Mollweide Plot
font = {'family': 'serif',
        'color':  'Black',
        'weight': 'normal',
        'size': 16,
        }
fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111, projection='mollweide')
Lon,Lat = np.meshgrid(lon,lat)
if scalesetting==0:
    minprob = arr1.min()   
    maxprob = arr1.max()
imm = ax.pcolor(Lon,Lat,arr1,norm=colors.Normalize(minprob,maxprob), cmap=plt.cm.jet, shading='auto')
cb = fig.colorbar(imm, extend='max', orientation='horizontal', shrink=0.67)
cb.set_label('P$_{a\gamma}$', fontdict=font)
ax.axis('off')
fig.savefig("skymap.png", format='png',dpi=300)
fig.savefig("skymap.pdf", format='pdf',dpi=100)
print('Process completed successfully!')
#Ref https://matplotlib.org/stable/api/figure_api.html#matplotlib.figure.Figure.colorbar
