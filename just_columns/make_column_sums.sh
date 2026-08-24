# Ale... nic... 

# Write header to column_sums.txt
echo -e "filename\t" > column_sums.txt
for n in 1, 2, 4, 8, 16, 32; do
    echo -e "top_${n}_sum\t" > column_sums.txt
done
echo -e "total_sum\t

for file in just_columns/*; do
    # Skip the output file itself if it's in the same directory
    [ "$file" = "just_columns/column_sums.txt" ] && continue
    [ "$file" = "just_columns/make_column_sums.sh" ] && continue


    
    for n in 1, 2, 4, 8, 16, 32; do
        # Calculate both sums in a single pass using awk
        awk -v cutoff="$n" '
            {
            row_sum = $1 + $2
            total_sum += row_sum
            if (NR <= cutoff) {
                top_n_sum += row_sum
            }
            }
            END {
            filename = FILENAME
            # Clean path to keep just the filename (e.g., ebola_n160_k31.tsv)
            sub(/^.*\//, "", filename)
            
            prop = (total_sum > 0) ? (top_n_sum / total_sum) : 0
            
            printf "%s\t%d\t%d\t%.4f\n", filename, top_n_sum, total_sum, prop
            }
        ' "$file" >> just_columns/column_sums.txt
    done
done