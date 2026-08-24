#!/bin/bash

# # Gemini

# # Define paired folders: (Core directory -> Result directory)
# declare -A folder_pairs=(
#     ["O_cores"]="original_results"
#     ["F_cores"]="fingerprint_results"
# )

# # Loop through each folder pair
# for core_dir in "${!folder_pairs[@]}"; do
#     result_dir="${folder_pairs[$core_dir]}"

#     echo "Processing pair: $core_dir -> $result_dir"

#     # Loop over all .fasta files in the core directory
#     for fasta_path in "$core_dir"/*.fasta; do
#         # Handle case where no files match the glob
#         [ -e "$fasta_path" ] || continue

#         # Extract base filename without directory or extension
#         filename=$(basename "$fasta_path" .fasta)
#         tsv_path="$result_dir/$filename.tsv"

#         # Check if corresponding .tsv file exists
#         if [ -f "$tsv_path" ]; then
#             # Count exact rows in the fasta file
#             row_count=$(wc -l < "$fasta_path")

#             # Simply append the count to the end of the .tsv file on a new line
#             echo "$row_count" >> "$tsv_path"

#             echo "  Appended core count ($row_count) to end of $filename.tsv"
#         else
#             echo "  Warning: Corresponding file not found: $tsv_path"
#         fi
#     done
# done

# echo "Done!"