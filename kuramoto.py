"""Goal here is to solve SDE in Uriu et al. (2020)."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# from process_tracks import get_data, get_neighbours

def kuramoto_coherence(array):
    """
    Calculates the kuramoto coherence order for an array of oscillators.

    Input: n x 1 array of phase values (floats).

    Output: float
    """
    z = 0
    for i in range(0, len(array)):
        z += np.exp(1j * array[i])
    z = z / len(array)
    return np.sqrt(z.real**2 + z.imag**2)

def deterministic(t, X, adjacency_array, omega=3, K=1):
    # omega = 10 # TODO Re-adjust this to random omega for each cell
    #  K = 1 # TODO Coupling constant, re-adjust this in future, or at least allow it to be adjusted
    """
    Calculates the first part of the Kuramoto SDE, for a set of oscillators.
    
    Input: 
    t (float) - time point.

    X (array) - array of phase values in the system.

    adjacency_array (array) - nxn array with 0 and 1 indicating if two oscillators are adjacent.
    Note that the oscillator is not considered adjacent to itself (entry neighbours[i,i]=0).

    omega (float) - intrinsic frequency of the cell. 

    K (float) - kuramoto coupling constant. 
    
    Output: n x 1 array of floats.
    """
    # Initialise the output array.
    output = np.zeros(len(X), dtype=float)
    # Assume that the order of X matches the order of the neighbouring array, and that they are the same size.
    assert len(X) == np.shape(adjacency_array)[0] and len(X) == np.shape(adjacency_array)[1], "Adjacency array should be square with edges of length equal to X."
    
    for index in range(0, len(X)):
        neighbours = adjacency_array[index]
        value = 0
        # Get the number of neighbours for this oscillator, and the indices of the oscillators in X
        number_of_neighbours = np.sum(neighbours)
        neighbour_indices = np.nonzero(neighbours)
        for array in neighbour_indices:
            for i in list(array):
                if i == index or i == 'NaN':
                    continue
                elif i != 'NaN':
                    value += np.sin(X[i]-X[index]) #IMPORTANT - index must not be in neighbour_indices
                    value = value * K / number_of_neighbours
                    value += omega
                    # Update array
                    output[i] = value 
                else:
                    raise ValueError(f"Index {i} in neighbour_indices array is not one of 0, 1, or 'NaN'.")
    return output 

def stochastic(time_step, num_vars, phase_noise=1.0): 
    """
    Returns a value for the stochastic term of the Kuramoto SDE.

    Input: 
    time_step (float) - dt
    num_vars (int) - number of variables in the set of SDEs.

    Output: array of floats
    """
    output = np.zeros(num_vars, dtype=float)
    for i in range(0, num_vars):
        delta_wt = np.random.normal(loc=0.0, scale=np.sqrt(time_step))
        output[i] = np.sqrt(2*phase_noise) * delta_wt
    return output

def euler_maruyama(time, X0, deterministic_part, stochastic_part, adjacency_array):
    """Numerically solves an SDE split into deterministic and stochastic components.
    
    Input:

    time (array) - array of time values for integration.
    X0 (array) - initial value for the integration.
    deterministic_part(t, X, adjacency_array) - array-valued function for deterministic component of the SDE.
    stochastic_part(time_step) - float-valued function for the stochastic component of the SDE

    Output: number_time_points x (number_variables+1) array. First row is the time value at each point,
    and remaining n rows are value of that variable at that time point. 
    """
    X = X0
    solution = np.zeros((len(X)+1, len(time)), dtype=float)
    for t in range(1, len(time)):
        solution[0, t] = time[t]
        dt = time[t] - time[t-1]
        #TODO handle errors due to divergence.
        X += dt * deterministic_part(time[t], X, adjacency_array) + stochastic_part(dt, len(X))
        for i in range(1, len(X)+1):
            solution[i, t] = X[i-1]
    return solution

# trial kuramoto by giving random phase
initial_condition = np.random.rand(10) # 10 oscillators
# give the time array
time = np.linspace(0, 100, 1000)
# adjacency_array (set all to be adjacent)
adjacency_array = np.ones((10,10))
for i in range(0, 10):
    adjacency_array[i,i] = 0 # so that not coupled to itself
# Integrate
solution = euler_maruyama(time, initial_condition, deterministic, stochastic, adjacency_array)
# plot the results
time_axis = solution[0]
# for i in range(0, len(initial_condition)):
#     y_axis = solution[i+1]
#     plt.plot(time_axis, y_axis)
# plot the coherence
coherence = np.zeros(len(time_axis))
for t in range(0, len(time_axis)):  
    time_slice = solution[1:,t]
    coherence[t] = kuramoto_coherence(time_slice)
plt.plot(time_axis, coherence)
plt.xlabel("Time")
plt.ylabel("Coherence of oscillators")
plt.show()


# Now we want to establish the initial phases of the simulation. 



# Now simulate the system 
