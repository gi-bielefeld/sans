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
typhimurium_path="/home/azetocha/typhimurium/fa"                     
ebola_path="/mnt/sans/ebola/fa/list.txt"    

# n_genomes:
a=$(wc -l < "$drosophila_path")
b=$(wc -l < "$ebola_path")

# ==============================================================================
# ORIGINAL SANS     BUILDING all SANS versions for all maxN: DONE
# ==============================================================================

# cd $sans_original

# for n in 1000 "$a" "$b"; do     	# omit 2000 3000 4000 5000, because that would be too much...
#     # Update maxN in the original makefile
#     sed -i "0,/DmaxN=/s/DmaxN=[0-9]*/DmaxN=${n}/" makefile
        
#     # FORCE CLEAN 
#     make -f "makefile" clean

#     # Build
#     make -f "makefile"

#     # Save the build of SANS under a unique name (using absolute paths):
#     mv "./SANS" "./SANS_n${n}"

#     echo "Build of the original SANS, n = ${n} was successful."
# done

# # Put maxN back as it was
# sed -i "0,/DmaxN=/s/DmaxN=[0-9]*/DmaxN=400/" makefile

# ==============================================================================
# FINGERPRINT SANS BUILD 
# ==============================================================================

cd $sans_fingerprint

# Typhimurium, 1000
# Set Fingerprint length to 64:
sed -i "0,/DFL=/s/DFL=[0-9]*/DFL=64/" makefile_F

for n in 1000; do        # omitting 2000 3000 4000 5000 
    # Update maxN
    sed -i "0,/DmaxN=/s/DmaxN=[0-9]*/DmaxN=${n}/" makefile_F
        
    # FORCE CLEAN the fingerprint directory
    make -f "makefile_F" clean

    # Build in the fingerprint directory
    make -f "makefile_F" 

    # Save the build of SANS under a unique name (using absolute paths):
    mv "./SANS" "./SANS_n${n}_FL64" 

    echo "Build of the fingerprint SANS, n = ${n}, FL = 64 was successful."
done


for n in "$a" "$b"; do 
    # Update maxN
    sed -i "0,/DmaxN=/s/DmaxN=[0-9]*/DmaxN=${n}/" makefile_F

    for FL in 16 32 64 128; do
        # Update Fingerprint length:
        sed -i "0,/DFL=/s/DFL=[0-9]*/DFL=${FL}/" makefile_F

        # FORCE CLEAN the fingerprint directory
        make -f "makefile_F" clean

        # Build in the fingerprint directory
        make -f "makefile_F" 

        # Save the build of SANS under a unique name (using absolute paths):
        mv "./SANS" "./SANS_n${n}_FL${FL}"

        echo "Build of the fingerprint SANS, n = ${n}, FL = ${FL} was successful."
    done
done 