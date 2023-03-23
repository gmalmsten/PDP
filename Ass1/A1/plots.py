# Read files and plot the data

import numpy as np
import matplotlib.pyplot as plt

# Read the data from the file
# Read from file in subfolder
input = np.loadtxt('input_file.txt', dtype=float, delimiter=' ')
output = np.loadtxt('out.txt', dtype=float, delimiter=' ')

input = input[1:]

fig, ax = plt.subplots()
ax.plot(input, 'b-', label='Input')
ax.plot(output, 'r-', label='Output')
ax.set_xlabel('Time')
ax.set_ylabel('Amplitude')
ax.legend()
plt.show()
