import sys
import os
import subprocess
import numpy as np
import pyphi
from phi_compute_function import *

# channel_sets should be specified in file 'networks_2ch', which holds channel_sets:
#	1,2
#	1,3
#	...
#	100,101
#	etc.

tpm_dir = sys.argv[1] # TPM directory (containing folder with TPMs, named same as output directory)
data_dir = sys.argv[2] # Source directory holding TPMs directory
out_dir = sys.argv[3] # Output directory
filename = sys.argv[4] # File holding TPMs

# Setup ############################################################################

pyphi.config.LOG_FILE = 'logs_pyphi/' + filename.split('.')[0] + '.log' # Log to file specific for this script

# Source directory
source_dir = tpm_dir + data_dir

# Actual causation (comment out for default IIT 3.0)
#pyphi.config.PARTITION_TYPE = 'ALL'
#pyphi.config.PICK_SMALLEST_PURVIEW = True

print(pyphi.config.PARTITION_TYPE + ' ' + str(pyphi.config.PICK_SMALLEST_PURVIEW), flush=True)

# Load data ############################################################################

print('loading data')

loaded = sio.loadmat(tpm_dir + data_dir + filename)
tpms = loaded['data']['tpms'][0][0] # Should be a matrix (states x nodes x networks)

print("loaded")

if len(tpms.shape) == 3:
	phis = np.zeros((tpms.shape[0], tpms.shape[2])) # states x networks
	state_sias = np.empty((tpms.shape[0], tpms.shape[2]), dtype=object)
	# Loop through
	for network in range(tpms.shape[2]):
		state_phis, state_sias = phi_compute(tpms[:, :, network])
		phis[:, network] = state_phis
		sias[:, network] = state_sias
else:
	phis = np.zeros((tpms.shape[0], 1))
	sias = np.empty((tpms.shape[0], 1), dtype=object)
	state_phis, state_sias = phi_compute(tpms[:, :])
	phis[:, 0] = state_phis
	sias[:, 0] = state_sias

# Save
sio.savemat(out_dir + filename, {'phis': phis, 'sias': sias}, do_compression=True, long_field_names=True)
print('saved ' + out_dir + filename, flush=True)
