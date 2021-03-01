import numpy as np
import pandas as pd


# As an alternative, write a function to extract the trajectory data to a numpy array
def get_data(path_to_data_file):
    """
    Extracts the .csv raw data to an mxnx5 array where each entry is for each cell at a given time.

    Input: .csv file with headers 'TrackID', 'Time', 'Position X Reference Frame', 
    'Position Y Reference Frame', 'Position Z Reference Frame'. 

    Output: num_cells x num_time_points x 5 numpy array.
    """
    print(f"\nReading raw data for file '{path_to_data_file}''...")
    data = pd.read_csv(path_to_data_file)
    # Now get the individual time points that are in the file
    times = list(set(data['Time']))
    # Get the cell identities
    cells = list(set(data['TrackID'])) # need to make it into a list, as sets don't support indexing
    # Now initialise the output array
    print("\nInitialising output array...")
    output = np.zeros((len(cells), len(times), 5))
    
    # Now we want to loop through the data and fill up the output array
    # Get length of each list for printing feedback.
    len_cells, len_times = len(cells), len(times)
    # Loop:
    for c_index, cell in enumerate(cells):
        # Subset by the cell's identity
        cell_data = data[data['TrackID'] == cell]
        for t_index, time in enumerate(times):
            # Give us live progress:
            print(f"\rWriting to output: cell {c_index}/{len_cells}, time {t_index}/{len_times}.", end="", flush=True)
            # Then for each specific time, get the coordinates, and write to the output array
            snapshot = cell_data[cell_data['Time'] == time]
            # need to check for the times where we have no data
            if snapshot.empty:
                continue
            else:
                x, y, z = snapshot[["Position X Reference Frame", 
                            "Position Y Reference Frame",
                            "Position Z Reference Frame"]].to_numpy(dtype=float)[0] # [0] important here - DON'T lose it
                # Get the coordinates for the output array
                col = times.index(time)
                row = cells.index(cell)
                # Now write to it
                output[row, col] = [float(time), x, y, z, cell]
    # Save to prevent calculating again
    output_file_path = path_to_data_file.replace(".csv", "") + "_as_trajectory.npy"
    np.save(output_file_path, output)
    print("\nData loaded.")
    return output


def get_neighbours(trajectory, file_name, r=10):
    """
    Returns an adjacency matrix with the neighbours for each cell in the simulation. 

    Input: 
    trajectory (array) - m x n x 5 array with cell positions over time. m = no. cells, n = no. time points
    file_name - name of the raw data file input 'trajectory' was extracted from.
    r (float) - radius in um in which to search for neighbouring nuclei.

    Output: t x n x n array with entries 0 and 1. Saved to file to minimise memory burden.

    Set default r as 8 in accordance with:
    Samsa LA, Fleming N, Magness S, Qian L, Liu J. 
    Isolation and Characterization of Single Cells from Zebrafish Embryos. 
    J Vis Exp. 2016;(109):53877. Published 2016 Mar 12. doi:10.3791/53877

    """
    # TODO maybe change r to be a function as input, to allow changes in r over time/space

    # first transpose the trajectory array so as to easily access individual time points
    print("\nGetting neighbours...")
    # Now get the number of cells
    num_cells = trajectory.shape[0]
    # And the number of time points in the data
    num_t = trajectory.shape[1]
    # now initialise the output matrix
    print("\nInitialising output...")
    output = np.zeros((num_t, num_cells, num_cells), dtype=int)
    # Now loop for each time slice and each cell, then check if each cell neighbours it
    # Now fill it up
    print("\nGenerating array...")
    for t in range(0, num_t):
        time_slice = trajectory[:,t,:]
        for i in range(0, num_cells):
            # Print for communication:
            print(f"\rGetting neighbours for cell {i}/{num_cells}, at time {t}/{num_t}.", end="")
            # get the cell's position
            xi, yi, zi = time_slice[i, 1:4]
            # Now calculate which cells are proximal
            for j in range(0, num_cells):
                # check if it's the same cell as above, in which case move on and set to 0
                if i == j:
                    continue
                # now check it's not a NaN value
                elif time_slice[j, 1] == 'NaN':
                    continue
                else:
                    # get the position of each nucleus and calculate if it's within r of the ith nucleus:
                    xj, yj, zj = time_slice[j, 1:4]
                    d = np.sqrt((xj - xi)**2 + (yj - yi)**2 + (zj - zi)**2) # euclidean distance
                    if d <= r:
                        # In this case, cells are adjacent
                        output[t,i,j] = 1
    # As the above is far too costly to execute repeatedly, save into a .npy file.
    output_name = "Adjacency_arrays/" + f"{file_name}_adjacency_array.npy"
    np.save(output_name, output)
    print(f"\nOutput array saved to file '{output_name}'.'")
    return output

# Now generate adjacency matrices for each track:
target_file = "M4_TrackingOnly.csv"
trajectory = get_data(target_file)
get_neighbours(trajectory, target_file)

# Now we want to interpolate the trajectories
"""This may now be not needed"""
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
    times = list(set(data['Time']))

    # Initialise the output array
    # time in rows, cell identity in columns, each entry a length 5 array with t, x, y, z coords and identifier
    n_cols = len(times) + n * (len(times)-1) 
    cells = list(set(data['TrackID'])) # need to make it into a list, as sets don't support indexing
    n_rows = len(cells) # get all the cell identifiers
    output = np.zeros((n_rows, n_cols, 5))

    # Now interpolate times, and create a large time array
    interpold_time = np.linspace(times[0], times[1], n, endpoint=False) # exclude the last element from the array
    for tx in range(1, len(times[2:])-1):
        interpold_time = np.append(interpold_time, np.linspace(times[tx], times[tx+1], n, endpoint=False))
    # Now need to add the final time value in the array
    interpold_time = np.append(interpold_time, [times[len(times)-1]])

    ### Now fill output up for each time point ### 

    # First enter in the values for track trajectories. These are needed to calculate interpolated coordinates.
    for cell in cells:
        # Subset by the cell's identity
        cell_data = data[data['TrackID'] == cell]
        for time in times:
            # Then for each specific time, get the coordinates, and write to the output array
            snapshot = cell_data[cell_data['Time'] == time]
            # need to check for the times where we have no data
            if snapshot.empty:
                continue
            else:
                x, y, z = snapshot[["Position X Reference Frame", 
                            "Position Y Reference Frame",
                            "Position Z Reference Frame"]].to_numpy(dtype=float)[0] # [0] important here - DON'T lose it
                # Get the coordinates for the output array
                col = np.where(interpold_time == float(time))
                row = cells.index(cell)
                # Now write to it
                output[row, col] = [float(time), x, y, z, cell]

    # Now for all other times, interpolate
    for c in range(0, n_rows):
        # Interpolate between the initial condition and the next non-zero coordinates
        # Find the indices of non-zero elements
        # Line below returns two arrays - one with coordinate along output[c], and the other coordinates in each entry
        # Only want the first one
        non_zeros = list(np.nonzero(output[c])[0])
        # Find the col (aka time) coordinate of each non-zero point
        # define coordinates vector as a set, so as to allow storage of many different time coordinates
        coords = set()  
        for non_zero in non_zeros:
            coords.add(non_zero) # column coordinate
        coords = list(coords) # to allow indexing
        # Now get each successive pair of coordinates
        for j in range(0, len(coords)-1):
            # matrix coordinates
            lower = coords[j]
            upper = coords[j+1]
            # time values
            t_lower = interpold_time[lower]
            t_upper = interpold_time[upper]
            # cartesian coordinates
            xl, yl, zl = output[c, lower][1:4] # TODO check for the NaN case after this
            xu, yu, zu = output[c, upper][1:4]
            X_l = [xl, yl, zl]
            X_u = [xu, yu, zu]
            # Now find out the time points we need to interpolate for, between t_lower and t_upper
            interpolations = interpold_time[lower+1:upper] # omit element with index 'lower', as it's specified by raw data
            # now for each time point between t_lower and t_upper, calculate the new x,y,z coordinates
            for t in range(0, len(interpolations)):
                time = interpolations[t]
                # Need a for loop to add vectors together to get the new coordinates
                new_coordinates = []
                for i in range(0, len(X_l)): # NOTE - assumes that len(Xl) == len(Xu) TODO - write a check before this
                    new_coordinate = X_l[i] + (X_u[i] - X_l[i]) * time / t_upper # equation of a straight line in 3D
                    new_coordinates.append(new_coordinate)
                x, y, z = new_coordinates
                output[c, t+lower+1] = [time, x, y, z, cells[c]]
    # Now change all 0 values to NaN, as these cells do not exist at this time point
    for i in range(0, n_rows):
        for j in range(0, n_cols):
            if output[i,j, 0]==0 and output[i,j, 1]==0 and output[i,j,2]==0 and output[i,j,3]==0 and output[i,j,4]==0:
                output[i,j] = np.array(['NaN', 'NaN', 'NaN', 'NaN', 'NaN'])
    # As the above takes forever to execute, save the numpy array for future use
    np.save("interpolated_tracks.npy", output)
    return output



# Now we want to get the cells in our system and define a neighbours array for each
# Will do this by an nxn matrix for the n cells in the system at time t

def get_neighbours_array(t, path_to_data_file, r=10):
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
