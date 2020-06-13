
import os
import sys
import psutil
import numpy as np
import scipy.io as sio
import pyphi
import math
import itertools

#https://pyphi.readthedocs.io/en/stable/configuration.html?highlight=after_computing_sia#pyphi.conf.PyphiConfig.CLEAR_SUBSYSTEM_CACHES_AFTER_COMPUTING_SIA
#pyphi.config.MAXIMUM_CACHE_MEMORY_PERCENTAGE = 50
# pyphi.config.CACHE_SIAS = False
# pyphi.config.CACHE_REPERTOIRES = False
# pyphi.config.CACHE_POTENTIAL_PURVIEWS = False
# pyphi.config.CLEAR_SUBSYSTEM_CACHES_AFTER_COMPUTING_SIA = False
pyphi.config.PROGRESS_BARS = False

# LEAVE AS FALSE - setting as true gives big slow down on MASSIVE
# Use TRUE if using multiple CPUs per task (SLURM), for searching many partitions
pyphi.config.PARALLEL_CUT_EVALUATION = False

def phi_compute(tpm):
	# Computes phi
	# Inputs:
	#	tpm = state-by-node matrix (states x nodes)
	# Outputs:
	
	# Build the network and subsystem
	# We are assuming full connection
	network = pyphi.Network(tpm)
	print("Network built", flush=True)

	#########################################################################################
	# Remember that the data is in the form a matrix
	# Matrix dimensions: sample(2250) x channel(15)

	# Determine number of system states
	n_states = tpm.shape[0]
	nChannels = int(math.log(n_states, 2)) # Assumes binary values

	# Initialise results storage structures
	state_sias = np.empty((n_states), dtype=object)
	state_phis = np.zeros((n_states))

	# sys.exit()

	# Calculate all possible phi values (number of phi values is limited by the number of possible states)
	for state_index in range(0, n_states):
		#print('State ' + str(state_index))
		# Figure out the state
		state = pyphi.convert.le_index2state(state_index, nChannels)
		
		# As the network is already limited to the channel set, the subsystem would have the same nodes as the full network
		subsystem = pyphi.Subsystem(network, state)
		
		#sys.exit()
		
		# Compute phi values for all partitions
		sia = pyphi.compute.sia(subsystem)
		
		#sys.exit()
		
		# Store phi and associated MIP
		state_phis[state_index] = sia.phi
		
		# MATLAB friendly storage format (python saves json as nested dict)
		# Store big_mip
		state_sias[state_index] = pyphi.jsonify.jsonify(sia)
		
		print('State ' + str(state_index) + ' Phi=' + str(sia.phi), flush=True)

	# Return ###########################################################################

	return state_phis, state_sias

########################################################################################################
########################################################################################################