#!/bin/bash

params=("$@")
# Name of the script to be generated
OUTPUT_SCRIPT1="spectrum.sh"
# Content of the generated script
cat <<EOL > "$OUTPUT_SCRIPT1"
#!/bin/bash
#SBATCH --partition=intelsr_devel
#SBATCH --time=00:50:00
#SBATCH --output=spectrumslurm-%j.out 
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=96

# Define the number of cores (or processes) to use
NUM_CORES=96

# Define the total number of tasks
TOTAL_TASKS=96

# Calculate the chunk size for each process
CHUNK_SIZE=\$((TOTAL_TASKS / NUM_CORES))

# Run each chunk in a separate process
for ((i=0; i<NUM_CORES; i++)); do
    START_IDX=\$((i * CHUNK_SIZE))
    if [ \$i -eq \$((NUM_CORES - 1)) ]; then
        END_IDX=\$TOTAL_TASKS  # Last chunk goes to the end
    else
        END_IDX=\$((START_IDX + CHUNK_SIZE))
    fi

    # Run the Python script for this chunk
    python3 spectrum.py ${params[@]} \$START_IDX \$END_IDX &
done

# Wait for all processes to finish
wait

# Combine all output files
python3 - <<EOF
import numpy as np
import os
import re

directory = "10/"

file_list = [f for f in os.listdir(directory) if f.endswith('.npy') and f.startswith('spectrum_')]

# Extract the numeric part and sort the files
file_list_sorted = sorted(file_list, key=lambda x: int(re.search(r'\d+', x).group()))

# Load and combine the arrays
combined_array = []
for file_name in file_list_sorted:
    file_path = os.path.join(directory, file_name)
    data = np.load(file_path)
    combined_array.append(data)

# Concatenate all arrays
combined_array = np.concatenate(combined_array, axis=1)

# Save the combined array to a new .npy file
output_file = os.path.join(directory, "spectrum.npy")
np.save(output_file, combined_array)

print(f"Combined data saved to {output_file}")
EOF

rm 10/spectrum_*
echo "All computations complete!"


TOTAL_TASKS=30
SECOND_NUM_CORES=30

# Calculate the chunk size for each process
CHUNK_SIZE=\$((TOTAL_TASKS / SECOND_NUM_CORES))

# Run each chunk in a separate process
for ((i=0; i<SECOND_NUM_CORES; i++)); do
    START_IDX=\$((i * CHUNK_SIZE))
    if [ \$i -eq \$((SECOND_NUM_CORES - 1)) ]; then
        END_IDX=\$TOTAL_TASKS  # Last chunk goes to the end
    else
        END_IDX=\$((START_IDX + CHUNK_SIZE))
    fi

    # Run the Python script for this chunk
    python3 Y.py ${params[@]} \$START_IDX \$END_IDX &
done

# Wait for all processes to finish
wait

# Combine all output files
python3 - <<EOF
import numpy as np
import os
import re

directory = "10/"

file_list = [f for f in os.listdir(directory) if f.endswith('.npy') and f.startswith('Y_')]

# Extract the numeric part and sort the files
file_list_sorted = sorted(file_list, key=lambda x: int(re.search(r'\d+', x).group()))

# Load and combine the arrays
combined_array = []
for file_name in file_list_sorted:
    file_path = os.path.join(directory, file_name)
    data = np.load(file_path)
    combined_array.append(data)

# Concatenate all arrays
combined_array = np.concatenate(combined_array, axis=0)

# Save the combined array to a new .npy file
output_file = os.path.join(directory, "Y.npy")
np.save(output_file, combined_array)

print(f"Combined data saved to {output_file}")
EOF

rm 10/Y_*
echo "All computations complete!"
EOL


chmod +x "$OUTPUT_SCRIPT1"

sbatch "$OUTPUT_SCRIPT1"
