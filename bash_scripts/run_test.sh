#!/bin/bash
# Exit immediately if any command fails (highly recommended for testing scripts)
set -e

# folder paths:
sans_original="/home/azetocha/SANS-master"
sans_fingerprint="/mnt/sans/SANS-fingerprinting2"
CURRENT_DIR=$(pwd)

mkdir -p original_info fingerprint_info                 
mkdir -p original_results fingerprint_results   
mkdir -p F_cores O_cores  

# data paths:
drosophila_path="/mnt/sans/drosophila/fa/wg_all.txt"            
typhimurium_path="/home/azetocha/typhimurium" 
ebola_path="/mnt/sans/ebola/fa/list.txt"                        

# ==============================================================================
# ORIGINAL SANS TESTING - DONE
# ==============================================================================

# cd $sans_original

# for data_set in "$ebola_path" "$drosophila_path"; do

#     n=$(wc -l < "$data_set")  
#     name=$(echo "$data_set" | cut -d/ -f4)
    
#     for k in 11 21 31; do

#         echo "========================================="
#         echo "Processing dataset: $name, k = ${k}, n = ${n}"
#         echo "========================================="  

#         # Define clean paths for outputs
#         output_file="original_results/${name}_n${n}_k${k}.tsv"
#         info_file="original_info/${name}_n${n}_k${k}.txt"

#         # Run and time original SANS											 	    Just one '>' to keep just last info!
#         { \time -v "$sans_original/SANS_n${n}" -i "$data_set" -r "O_cores/${name}_n${n}_k${k}.fasta"  -R "$output_file" -v -k "${k}" -T 32; } > "$info_file" 2>&1

#     done
# done 

# # Typhimurium
# for n in 1000; do       # omit  2000 3000 4000 5000

#    data_set="$typhimurium_path/list_${n}.txt"

#     for k in 11 21 31; do      # k = 11 and 21 was done already

#         echo "========================================="
#         echo "Processing dataset: Typhimurium (sub_${n}), k = ${k}"
#         echo "========================================="    

#         # Define clean paths for outputs
#         output_file="original_results/typhimurium_n${n}_k${k}.tsv"
#         info_file="original_info/typhimurium_n${n}_k${k}.txt"
        
#         # Run and time original SANS
#         { \time -v "$sans_original/SANS_n${n}" -i  "$data_set" -r "O_cores/typhimurium_n1000_k${k}.fasta" -R "$output_file" -v -k $k -T 32; } > "$info_file" 2>&1
    
#     done
# done

# ==============================================================================
# FINGERPRINT SANS TESTING
# ==============================================================================

cd $sans_fingerprint

for k in 11 21 31; do     

    for FL in 16 32 64 128; do 

        # Drosophila and Ebola
        for data_set in "$ebola_path" "$drosophila_path"; do

            # Recalculate 'n' specifically for these datasets to match the pre-built binary
            n=$(wc -l < "$data_set")  
            name=$(echo "$data_set" | cut -d/ -f4)

            echo "========================================="
            echo "Processing dataset: $name, k = ${k}, FL = ${FL}, n = ${n}"
            echo "========================================="
            
            output_file="fingerprint_results/${name}_n${n}_k${k}_FL${FL}.tsv"
            info_file="fingerprint_info/${name}_n${n}_k${k}_FL${FL}.txt"

            # Use the specific pre-built SANS
            { \time -v "$sans_fingerprint/SANS_n${n}_FL${FL}" -i "$data_set" -R "$output_file" -r "F_cores/${name}_n${n}_FL${FL}_k${k}.fasta" -v -k $k -T 32; } > "$info_file" 2>&1
        done
    done

    FL=64

    # Typhimurium sublists
    for n in 1000; do      # omit  2000 3000 4000 5000
        data_set="${typhimurium_path}/list_${n}.txt"

        echo "========================================="
        echo "Processing dataset: Typhimurium (sub_${n}), k = ${k}, FL = ${FL}"
        echo "========================================="
        
        output_file="fingerprint_results/typhimurium_n${n}_k${k}_FL${FL}.tsv"
        info_file="fingerprint_info/typhimurium_n${n}_k${k}_FL${FL}.txt"
        
        # Use the specific pre-built SANS
        { \time -v "$sans_fingerprint/SANS_n${n}_FL${FL}" -i "$data_set" -R "$output_file" -r "F_cores/typhimurium_n${n}_FL${FL}_k${k}.fasta" -v -k $k -T 32; } > "$info_file" 2>&1
    done
done
