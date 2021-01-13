"""Goal here is to solve SDE in Uriu et al. (2020), before translating to Julia."""

import numpy as np
import pandas as pd

def deterministic(t, X, index, neighbours):
    """
    Calculates the first part of the SDE, for a specific oscillator.
    
    Input: 
    t (float) - time point.
    X (array) - array of phase values in the system.
    index (int) - position of the oscillator in the array X.
    neighbours (array) - nx1 array with 0 and 1 indicating if two oscillators are adjacent.
    Note that the oscillator itself is within this array but is given as 0.

    Output: float.
    """
    output = 0
    omega = 10 # TODO Re-adjust this to random omega for each cell
    K = 1 # TODO Coupling constant, re-adjust this in future, or at least allow it to be adjusted
    # Get the number of neighbours for this oscillator, and the indices of the oscillators in X
    number_of_neighbours = np.sum(neighbours)
    neighbour_indices = np.nonzero(neighbours)
    # TODO put a check here that index is not in neighbour_indices
    for i in neighbour_indices:
        output += np.sin(X[i]-X[index]) #IMPORTANT - index must not be in neighbour_indices
    output = output * K / number_of_neighbours
    output += omega
    return output 

def stochastic(time_step, phase_noise=1.0): 
    """
    Returns a value for stochastic term.

    Input: time_step (float) - dt

    Output: float
    """
    delta_wt = np.random.normal(loc=0.0, scale=np.sqrt(time_step))
    output = np.sqrt(2*phase_noise) * delta_wt
    return output

# Now we want to interpolate the trajectories
def interpolate_tracks(path_to_data_file, n=100):
    """
    Interpolates the cell movement between timepoints.

    Inputs: 
    path_to_data_file (str) - path to .csv file containing cell tracking data.
    n (int) - number of time points generated between each two time points in the raw data.

    Output: array with cell positions at specified time
    """
    # Need the position of each cell at successive time points
    data = pd.read_csv(path_to_data_file)
    times = set(data['Time'])

    # Initialise the output array
    # time in rows, cell identity in columns, each entry a length 5 array with t, x, y, z coords and identifier
    n_cols = len(times) + n * (len(times)-1) 
    cells = list(set(data['TrackID'])) # need to make it into a list, as sets don't support indexing
    n_rows = len(cells) # get all the cell identifiers
    output = np.zeros((n_rows, n_cols, 5))

    # Now interpolate times, and create a large time array
    interpold_time = np.linspace(times[0], times[1], n, endpoint=False) # exclude the last element from the array
    for tx in range(0, len(times[2:])-1):
        interpold_time = np.append(interpold_time, np.linspace(times[tx], times[tx+1], n, endpoint=False))
    # Now need to add the final time value in the array
    interpold_time = np.append(interpold_time, [times[len(times)]])

    ### Now fill output up for each time point ### 

    # First enter in the values for track trajectories. These are needed to calculate interpolated coordinates.
    for cell in cells:
        # Subset by the cell's identity
        cell_data = data[data['TrackID'] == cell]
        for time in times:
            # Then for each specific time, get the coordinates, and write to the output array
            snapshot = cell_data[cell_data['Time'] == time]
            x, y, z = cell_data["Position X Reference Frame", 
                        "Position Y Reference Frame",
                        "Position Z Reference Frame"].to_numpy(dtype=float)
            # Get the coordinates for the output array
            col = np.where(interpold_time == float(time))
            row = cells.index(cell)
            # Now write to it
            output[row, col] = [float(time), x, y, z, cell]

    # Now for all other times, interpolate
    for c in range(0, n_rows):
        # Interpolate between the initial condition and the next non-zero coordinates
        # Find the indices of non-zero elements
        non_zeros = np.non_zero(output[c]) # This returns 2-tuples
        # Find the row, col coordinates of each non-zero point
        i = 0 
        coords = []
        while i <= len(non_zeros): # TODO first check that len(non_zeros) % 5 == 0
            coords.append(non_zeros[i][0])
            i += 5  # NOTE - we assume there are ALWAYS 5 elements in each non-zero entry.
        # Now get each successive pair of coordinates from coords, and interpolate
        for j in range(1, len(coords)+1):
            # matrix coordinates
            lower = coords[j-1]
            upper = coords[j]
            # time values
            t_lower = output[]


            









# Now we want to get the cells in our system and define a neighbours array for each
# Will do this by an nxn matrix for the n cells in the system at time t

def get_neighbours_array(t, path_to_data_file, r=20):
    """
    Returns an adjacency matrix with the neighbours for each cell in the simulation. 

    Input: 
    t (int) - time step at which to get the neighbours array.
    path_to_data_file (str) - path to .csv file where the tracking data is stored.
    r (float) - radius in um in which to search for neighbouring nuclei.

    Output: nxn array with entries 0 and 1. 
    """
    # TODO take time x cells x 3 array as input instead
    # TODO maybe allow plotting of the connectivity graph over time?
    # Load the dataframe
    data = pd.read_csv(path_to_data_file)
    # Subset the dataframe by the time value
    data = data[data["Time"] == t]
    # Now get the list of cell IDs
    cells = set(data["TrackID"])
    # NOTE - the ordering of the cells in the neighbours array is determined by the ordering
    # of the 'cells' set above. 
    # Define the output
    # Note that neighbours for each cell are stored in rows.
    output = np.zeros((len(cells), len(cells)))
    # Now loop over each cell and calculate the neighbours for each
    i = 0
    for celli in cells:
        # Get the cell's position
        position_data = data[data["TrackID"] == celli]
        # NOTE - all cell IDs should be in data, as we subset by time point
        # Extract the coordinates from a numpy array
        xi, yi, zi = position_data["Position X Reference Frame", 
                "Position Y Reference Frame",
                "Position Z Reference Frame"].to_numpy()
        # Now loop over all cells and check if they're adjacent
        j = 0
        for cellj in cells:
            if cellj == celli:
                # assume every cell is NOT a neighbour of itself
                continue 
            else:
                # Get the coordinates of that cell
                position_data = data[data["TrackID"] == cellj]
                xj, yj, zj = position_data["Position X Reference Frame", 
                    "Position Y Reference Frame",
                    "Position Z Reference Frame"].to_numpy()
                # Check if it's adjacent using euclidean metric
                d = np.sqrt((xj - xi)**2 + (yj - yi)**2 + (zj - zi)**2)
                if d <= r:
                    # In this case, cells are adjacent
                    output[i,j] = 1
            # update cell j index
            j += 1
        # update cell i index
        i += 1
    return output

# Now simulate the system 
