# Adrian Zetocha        23.6.2026
# 
# This script extracts performance information from info and split files.

import subprocess
import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
from math import ceil

cols = ["dataset", "version", "k", "n", "FL", "running time", "CPU percent", "running time adjusted", "max memory (kb)", 
        "kmers read", "singletons percentage", "kmers per color set",
        "number of color sets", "precision", "threads"]

species = {"typhimurium":[1000], "drosophila":[12], "ebola":[160]}

# Also it would be interesting to pinpoint the distribution of kmers in color sets:
# show that most are concentrated in a top few.

df = pd.DataFrame(columns=cols);

# Good practice: always ensure all the folders the program needs exist.
subprocess.run("mkdir -p result_discrepancies just_columns" , shell=True,)


def parse_time_string(time_str):
    "Outputs seconds."
    if time_str[0] == '(':
        time_str = time_str[1:]
    if time_str[-1] == ')':
        time_str = time_str[:-1]

    number, unit = time_str.split()
    multiple = 0
    if (unit == "ms"):
        multiple = 1 / 1000
    elif (unit == "sec"):
        multiple = 1
    elif (unit == "min"):
        multiple = 60
    else:
        raise(ValueError("Invalid unit in parse time string: ", time_str))
    
    return float(number) * multiple
    

def read_results(file_name):

    # This part counts the number of non-zero count color sets from the splits file: (just_columns)
    n_color_sets = 0
    try:
        result = subprocess.run("awk -F' ' '{if ($1 > 0) {n += 1}; if ($2 > 0) n += 1} END" \
                                " {print n}' " +  f"{file_name}", 
                                shell=True, capture_output=True, text=True, check=True)

        # Print the command's output
        # print("Command Output:")
        # print(result.stdout)
        n_color_sets = int(result.stdout)

    except subprocess.CalledProcessError as e:
        print(f"Error occurred while executing the command: {e}")
        print(f"Error Output: {e.stderr}")


    # Determine precision
    precision = 100
    if "FL" in file_name:
        # fingerprinting version - we need to compare
        # file_name is has form: fingerprint_results/species_n123_k123_FL123.tsv
        
        # Multi-threading tests:
        # works still!

        # keep just the file name:
        ff_without_folder = file_name[file_name.index('/')+1 : ]
        # FL is the last parameter - prune it (in case of multithreading tests it's the pre-last)
        original_file = "original_results/" + ff_without_folder[0 : ff_without_folder.index("_FL")] + ".tsv"
        # original_file = original_results/species_n123_k123.tsv
        diff_file = "result_discrepancies/" + ff_without_folder

        # already compared and saved comparison to diff file.        
        command = f"( diff --suppress-common-lines -y {original_file} {file_name} ) > {diff_file} 2>&1"
        command_result = subprocess.run(command, shell=True)
        print(f"return code of command: {command}:\n", command_result.returncode )
        print(f"diff file tail: \n", subprocess.run(f"tail -n 10 {diff_file}", shell = True, capture_output=True, text=True, check=True))

        # check the diff file:
        with open(diff_file, "r") as file:
            if (file.readline() == ""):
                precision = 100
                print("The results in: ", ff_without_folder, "are perfect.")
            else:
                precision = "Nan"
                print("The results in: ", ff_without_folder, "differ from the original.")
            
    return precision, n_color_sets


def read_info(info, results, species, n, k, FL = None):
    try:
        with open(info, 'r') as file:
            pass
    except FileNotFoundError:
        print("File: ", info, "was not found. Skipping.")
        return
    with open(info, 'r') as file:
        # ... kmers read"
        line = file.readline()
        while not(ord('0') < ord(line[0]) < ord('9')): 
            line = file.readline()
        kmers_info = line.split()
        assert(kmers_info[1] == "k-mers" and kmers_info[2] == "read.")
        kmers_read = int(kmers_info[0])
        singleton_percentage = int(kmers_info[5][0:-1])

        # Done! ( ... ms)"
        while line[0:6] != " Done!":
            line = file.readline()

        time_string = line[7:-1]  # without brackets
        exec_time_s = parse_time_string(time_string)

        very_next_line = file.readline()
        assert("Command being timed:" in very_next_line)
        words = very_next_line.split()
        thr = words[words.index("-T") + 1]
        if '"' == thr[-1]:
            thr = thr[0:-1]
        N_THREADS = int(thr)

        while "Percent of CPU this job got:" not in line:
            line = file.readline()

        # Added CPU_percent because we should account for it
        # 500% means an equivalent of 5 cores on 100% were used during this job.
        CPU_percent = int(line.split()[-1][:-1])

        while "Maximum resident set size (kbytes):" not in line:
            line = file.readline()

        max_memory_consumption = line.split()[-1]

        version = "sans original"
        if FL is not None:
            version = "fingerprint"

    kpcs = 0
    precision, n_color_sets = read_results(results)
    if n_color_sets == 0:
        print("Error: zero color sets!")
        print("Current File:  ", info, results)
    else: 
        kpcs = round(kmers_read / n_color_sets)

    # cols = ["dataset", "version", "k", "n", "FL", "running time", "max memory (kb)", 
    # "kmers read", "singletons percentage", "kmers per color set",
    # "number of color sets", "precision"]

    new_row = { "dataset" : species,
                "version" : version,
                "k" : k,
                "n" : n,
                "FL" : FL,
                "running time" : exec_time_s,
                "CPU percent" : CPU_percent,
                "running time adjusted" : exec_time_s * (CPU_percent) // 100,
                "max memory (kb)" : max_memory_consumption, 
                "kmers read" : kmers_read,
                "singletons percentage" : singleton_percentage,
                "kmers per color set" : kpcs,
                "number of color sets" : n_color_sets,
                "precision" : precision,
                "threads" : N_THREADS
                }

   # save to df
    df.loc[len(df)] = new_row


# subprocess.run("mkdir -p just_columns" ) 

def write_kmer_distribution(file_name):
    # Distribution of kmers in color sets.
    # Remove the color set labels (keep just the two number columns)
    keep_just_columns = "awk -F'\t' '{print $1, $2}'"
    save_name = f"just_columns/{file_name}"
    subprocess.run( f"{keep_just_columns} {file_name} > {save_name}")

# read data 
if (False):
    for s in species.keys():
        for k in (11, 21, 31):
            for n in species[s]:
                # ORIGINAL RESULTS
                res_orig =  f"original_results/{s}_n{n}_k{k}.tsv"
                info_orig = f"original_info/{s}_n{n}_k{k}.txt"

                read_info(info_orig, res_orig, s, n, k)
                # write_kmer_distribution(res_orig)         # done

                # FINGERPRINTING RESULTS 
                for FL in (16, 32, 64, 128):
                    res_F =  f"fingerprint_results/{s}_n{n}_k{k}_FL{FL}.tsv"
                    info_F = f"fingerprint_info/{s}_n{n}_k{k}_FL{FL}.txt"

                    print("\n = = = = = = = = Processing species:", s, "n =", n, "k =", k, "FL =", FL, " = = = = = = = =")
                    read_info(info_F, res_F, s, n, k, FL)

# Testing threads:
if (True):
    for s in species.keys():
        for k in (11, ):
            for n in species[s]:
                for T in (1, 2, 3, 5, 10, 16):
                    # # ORIGINAL RESULTS
                    # res_orig =  f"original_results/{s}_n{n}_k{k}.tsv"
                    # info_orig = f"original_info/{s}_n{n}_k{k}.txt"

                    # read_info(info_orig, res_orig, s, n, k)
                    # # write_kmer_distribution(res_orig)         # done

                    # FINGERPRINTING RESULTS - with threads and only FL 64
                    for FL in (64, ):
                        res_F =  f"thread_results/{s}_n{n}_k{k}_FL{FL}_T{T}.tsv"
                        info_F = f"thread_info/{s}_n{n}_k{k}_FL{FL}_T{T}.txt"

                        print("\n = = = = = = = = Processing species:", s, "n =", n, "k =", k, "FL =", FL, "threads =", T, "= = = = = = = =")
                        read_info(info_F, res_F, s, n, k, FL)

# read_info("vzor_vysledku.txt", "vysledok_classic_f.tsv", "virus", 13, 11)
df.to_csv('TABLE_multithreading.csv', index = False)