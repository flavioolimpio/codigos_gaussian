import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

# Define a function to calculate the step size
def calculate_step(y_min, y_max, num_ticks):
    # Calculate the range of the data
    y_range = y_max - y_min
    
    # Calculate the approximate step size
    step = y_range / (num_ticks - 1)
    
    # Round the step size to a nice value
    scale = 10 ** np.floor(np.log10(step))
    nice_step = scale * np.round(step / scale)
    
    return nice_step

# Set the directory containing the directories to be entered
#dir = "C:\\Users\\Flavio\\OneDrive\\Flavio\\Pos_Doc_UEG\\Colaboracoes\\Elisa\\Calculos\\Fukui\\"

dir = "C:\\Users\\flavi\\OneDrive\\Flavio\\Pos_Doc_UEG\\Colaboracoes\\Elisa\Calculos\\Fukui\\"

# Create a list of directory names
dirs = ["Amoxicillin", "Metronidazole", "Cefazolin", "Sulfamethoxazole", "Chloramphenicol"]

# Create a figure with a GridSpec layout
fig = plt.figure(figsize=(15, 10))
gs = GridSpec(nrows=2, ncols=6, figure=fig, wspace=0.6)

# Loop through the directory names
for i, d in enumerate(dirs):
    # Read the CDFT.txt file and remove the first 4 lines
    with open(dir + d + '\\' + "CDFT.txt", "r") as f:
        cdft_lines = f.readlines()[5:]
    
    # Read the N.wfn file and select the number in the penultimate column of the second line
    with open(dir + d + '\\' + "N.wfn", "r") as f:
        n_wfn_lines = f.readlines()
        num_lines = int(n_wfn_lines[1].split()[-2])
    
    # Use the selected number to determine the number of lines to read from CDFT.txt
    atoms = [line.split()[0].replace('(','_').replace(')', '') for line in cdft_lines[:num_lines]]
    f0 = [float(line.split()[-2]) for line in cdft_lines[:num_lines]]

    data = {'Atoms' : atoms, 'Fukui' : f0 }
    
    df = pd.DataFrame(data)
    df.sort_values(by='Fukui', ascending=False, na_position='first', inplace=True)
    
    # Determine the row and column index of the subplot
    row = i // 3
    col = (i % 3) * 2 + (i // 3) * 1
    
    # Create an Axes object for the subplot
    ax = fig.add_subplot(gs[row, col:col+2])
    
    # Create the bar plot in the corresponding subplot
    ax.bar(df['Atoms'][0:5], df['Fukui'][0:5], color ='blue', width = 0.4)
     
    ax.set_xlabel("Atoms in  " + d + " molecule")
    ax.set_ylabel("Fukui Parameter (f0)")

    # Calculate and set the step size of the yticks
    y_min, y_max = ax.get_ylim()
    num_ticks = 5
    step = calculate_step(y_min, y_max, num_ticks)
    
    ax.set_yticks(np.arange(y_min, y_max + step, step))
    
    # Create a new list of labels in the desired format
    labels = [f"{atom.split('_')[1]}({atom.split('_')[0]})" for atom in df['Atoms']]

    # Set the xtick labels to the new list of labels
    ax.set_xticklabels(labels)
    
    
# Display the figure
plt.show()
#plt.savefig("C:\\Users\\flavi\\OneDrive\\Flavio\\Pos_Doc_UEG\\Colaboracoes\\Elisa\\Figuras\\" + "Figure_Fukui.png", format='png', dpi=600)
