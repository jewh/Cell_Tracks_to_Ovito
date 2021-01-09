"""
This script visualises cell tracking data. More info to follow. 
Email james.hammond@merton.ox.ac.uk for further info. 

Brief idea:
'Cell' Class Object - which has position attributes, id, etc. 
Use this to access things about a given cell

Plot protocol - a function or similar which plots all the cells over time, and animates them. 
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ovito.io as ov

class Cell:
    # TODO maybe remove the time argument and allow us to call a track for a single cell?
    def __init__(self, identifier, path_to_data):
        self.id = identifier
        self.path_to_data = path_to_data

    def trajectory(self):
        """
        Gets the position of the cell over time.
        
        Returns a numpy array. 

        """
        if str(self.path_to_data).find(".xlsx") != -1: # TODO make this check only last 5 chars
            data = pd.read_xlsx(self.path_to_data)
        elif str(self.path_to_data).find(".csv") != -1:
            data = pd.read_csv(self.path_to_data)
        else:
            raise ValueError(f"'{self.path_to_data}' is not an .xlsx or .csv file.")
        # Now subset the data for the given id and time
        data = data.loc[data['TrackID'] == self.id] 
        # Now extract the coordinates of the cell's position at each time point
        # (Convert to numpy array to permit a for loop).
        data = data.to_numpy(na_value = 'NAN') 
        trajectory = np.zeros((len(data), 4), dtype=float)
        for i in range(0, len(data)):
            row = data[i]
            x, y, z, t = row[0], row[1], row[2], row[4]
            trajectory[i][0:4] = x, y, z, t
        return trajectory


    def position(self, time):
        """Reads the cell's position at the given time from the trajectory."""
        for row in self.trajectory():
            if row[3] == float(time):
                x, y, z = row[0:3]
                return x, y, z
        return 'NaN', 'NaN', 'NaN'

def plot_cells(time, path_to_data):
    # First get the cell IDs
    ids = set(pd.read_csv(path_to_data)['TrackID'])
    # Now for each of these ids get the position at that time
    positions = np.zeros((len(ids), 3))
    for i in ids:
        cell = Cell(i, path_to_data)
        position = cell.position(time)
        del cell # save memory
        positions[i][0:3] = position
    return positions

data_file = "M4_TrackingOnly.csv"

# unwrapping the track data into snapshots, then save these, then access using ovito
data = pd.read_csv(data_file)
times = set(data['Time'])
# save a snapshot for each time
for t in times:
    snapshot = data[data["Time"] == t] # subset by that time
    snapshot = snapshot[["Position X Reference Frame", 
                    "Position Y Reference Frame",
                    "Position Z Reference Frame", 
                    "TrackID"]] # subset this by relevant columns
    snapshot = snapshot.to_numpy()
    # save to .txt
    np.savetxt(f"Snapshots/snapshot_t{t}.txt", snapshot, fmt="%.6f", delimiter = "\t") 
    # TODO note the number of decoimal places seems to be different from the .csv file and I can't fix that
