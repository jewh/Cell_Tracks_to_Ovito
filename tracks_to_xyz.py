"""
Aim is to convert .csv cell tracking data to .xyz format so as to plot as an ovito animation.
"""

import numpy as np
import pandas as pd

def to_xyz(array, time):
    """
    Input: Array with four columns. Assumes the first three are cartesian coordinates. 
    Fourth is assumed to be some property of the cells at each coordinate - identity, 
    level of gene expression, etc. 

    Output: saves input array as an .xyz file
    """
    # Get the number of 'atoms'. This is the no. of particles in the .xyz file
    no_atoms = str(len(array))
    # Get the comment for the file (goes on line 2).
    comment = f"Snapshot for time = {time}."
    # initialise the output file
    file_name = f"Snapshots_genes/snapshot_t{time}.xyz"
    with open(file_name, "w") as f:
        # write the first two lines
        f.write(no_atoms + "\n" + comment + "\n")
        # TODO allow tracking of cells by identity 
        # perhaps let them take on an identity depending on simulation?
        # Now loop over each set of coords in the array and append to the .xyz file
        for row in array:
            x, y, z, identity = row
            f.write(str(x) + "\t" + str(y) + "\t" + str(z) + "\t" + str(identity) + "\n")
    
data_file = "gene_expression_M4_TrackingOnly.csv"

# Now convert to .xyz
# Subset by time 
data = pd.read_csv(data_file) # TODO alter script to accept .txt or .xlsx files
times = set(data['Time'])
# save a snapshot for each time
for t in times:
    snapshot = data[data["Time"] == t] # subset by that time
    snapshot = snapshot[["Position X Reference Frame", 
                    "Position Y Reference Frame",
                    "Position Z Reference Frame", 
                    "gene_expression"]] # subset this by relevant columns
    snapshot = snapshot.to_numpy()
    to_xyz(snapshot, t)



    