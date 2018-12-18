# Module that collects frequencies, amplitudes and phases from specified files
# to form one input file for the pso program

# It requires that pso program uses single file containing frequencies, amplitudes and phases in adequete columns:
# freq  amp  phase

# The last parameter specifies how many first amplitude samples should be truncated to 0

import sys

import numpy as np

input_files = [
    sys.argv[1],
    sys.argv[3],
    sys.argv[5]
]
input_channels = [
    int(sys.argv[2]),
    int(sys.argv[4]),
    int(sys.argv[6])
]
output_file = sys.argv[7]

smooth_samples = int(sys.argv[8])

print(input_files)
print(input_channels)
print(output_file)

data = []

for i in xrange(len(input_files)):
    input_file = input_files[i]
    column = input_channels[i]
    with open(input_file) as f:
        data.append([float(r.split()[column]) for r in f.readlines()])

# Smooth specified number of first samples
for i in xrange(smooth_samples):
    data[1][i] = 0.0

np.savetxt(
    output_file,
    np.transpose(data),
    delimiter='\t'
)
