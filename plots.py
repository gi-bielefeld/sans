# Adrian Zetocha        23.6.2026

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re


os.makedirs("plots", exist_ok=True)

df = pd.read_csv("TABLE.csv")
df.rename(columns={"dataset": "species",
                   "max memory (kb)": "max memory (GB)",
                   "running time": "running time (min)",
                   "running time adjusted" : "running time (min, adjusted)"}, inplace=True)


# Changed run time from seconds to minutes:
df['running time (min)'] /= 60
df['running time (min, adjusted)'] /= 60
df['max memory (GB)'] /= 10**6

# In numpy, slightly different logiacl operators are used!: | for or
df['version'] = np.where((df['precision'] == 100) | (df['precision'] == "100"),
                         df['version'], "incorrect results")

# print("Head of df - for multithreading:")
print("Head of TABLE.csv:")
print(df.head())
species = {"typhimurium": [1000, 2000], "drosophila": [12], "ebola": [160]}

# ===== ==== === For multithreading = = = == = = = = = = == = = 

# df_prev = pd.read_csv("TABLE.csv")
# df_prev.rename(columns={"dataset": "species",
#                    "max memory (kb)": "max memory (GB)",
#                    "running time": "running time (min)",
#                    "running time adjusted" : "running time (min, adjusted)"}, inplace=True)


# # Changed run time from seconds to minutes:
# df_prev['running time (min)'] /= 60
# df_prev['running time (min, adjusted)'] /= 60
# df_prev['max memory (GB)'] /= 10**6

# df_prev['version'] = np.where(df_prev['precision'] == '100',
#                          df_prev['version'], "incorrect results")


# # For plotting
# # Nezabudni pridat tri riadky pre T = 32 z TABLE.csv
# temp = df_prev.query("FL == 64 and k == 11")
# temp["threads"] = 32            # v celom stlpci bude rovnaka hodnota 32
# df = pd.concat([df, temp])           # pridaj tabulku k prvej

# # a prazdne riadky pre T = 4, 6-9, 11-15, 17-31.
# for t in (4, 6,7,8,9,11,12,13,14,15,17):
#     print("df has:", len(df.index), "rows.")
#     df.loc[len(df.index)].all() = pd.DataFrame([None for _ in range(df.shape[1] -1)] + [t] )

# for t in range(18, 33):
#     df.loc[len(df.index)].all() = pd.DataFrame([None for _ in range(df.shape[1] -1)] + [t] )


def format_sci(val):
    formatted = f"{val:.2e}"  # Produces e.g. '1.45e+08'
    base, exp = formatted.split('e+')
    return f"{base} x 10^{int(exp)}"


def plot_all_species_old(species):

    # Set clean visual style
    sns.set_theme(style="whitegrid")
    plt.tight_layout()
    plt.xticks((11,21,31))

    row = 0; col = 0; nrows = 2; ncols = 2

    # comparison: sans original - fingerprint:  "FL is None or FL == 64"
    # comparison between fingerprints:  "FL is not None"
    for counter, condition in enumerate(("FL.isna() or FL == 64", "not FL.isna()")):          # 'Is' nodes are not implemented
        for parameter in ("running time", "max memory (kb)"):

            fig, axes = plt.subplots(nrows = nrows, ncols = ncols, figsize=(16, 20), sharey=False, sharex=False)
            fig.suptitle(parameter)

            for name, counts in species.items():
                for n_genomes in counts:
                    species_df = df.query(f"({condition}) and (species == '{name}') and (n == {n_genomes})")
                    ax = axes[row][col]
                    # increment
                    col += 1
                    if col == ncols:
                        col = 0; row += 1
                    if row == nrows:
                        row = 0
                    
                    print(name, "df :")
                    print(species_df)
                    
                    if species_df.empty:
                        # ax.text(0.5, 0.5, f"No data for {species}", ha='center', va='center')
                        # ax.set_title(species.capitalize())
                        continue
                        
                    # Sort logically by k so the columns go left-to-right (e.g., k=11, k=21, k=31)
                    species_df = species_df.sort_values(by='k')
                
                    # 3. GENERATE THE SCATTER PLOT
                    # We use hue to split colors by "version" (Original vs Fingerprint)
                    sns.scatterplot(
                        data=species_df,
                        x ='k',
                        y=parameter,
                        hue='version',
                        style='version',   # Gives different markers (e.g., circle vs X) for extra clarity
                        s=150,             # Distinct, large data points
                        ax=ax,
                        palette={'sans original': 'green', 'fingerprint': 'blue'}, 
                        edgecolor='none',
                        alpha=0.85
                    )

                    # 4. ANNOTATE DATA POINTS
                    for _, r0w in species_df.iterrows():
                        
                        # Format annotation: display precision score or FL if applicable
                        if r0w['version'] == 'fingerprint':
                            label_text = f"FL{int(r0w['FL'])}"
                        else:
                            label_text = "sans original"

                        # drop text labels right next to data points
                        ax.text(
                            x = r0w['k'] + 2, 
                            y = r0w[f"{parameter}"]-5,
                            s=label_text, 
                            fontsize=13, 
                            ha='center', 
                            va='bottom',
                            weight='semibold'
                        )

                    # 5. Annotate x ticks
                    tick_positions = (11,21,31)
                    tick_labels = []

                    # let's take the average of kmers read for each k
                    for k in (11,21,31):
                        column = species_df.query(f"k == {k}")['kmers read']
                        avg = column.sum() / len(column)
                        s = f"{k}\n" + str(format_sci(avg))
                        tick_labels.append(s)

                    # Set tick positions first
                    ax.set_xticks(tick_positions)

                    # Set labels second
                    ax.set_xticklabels(tick_labels, rotation=0, fontsize=10)

                    # 5. REFINING AXES AND TITLES
                    ax.set_title(f" {name}, n = {n_genomes}", fontsize=14, pad=10, weight='bold')
                    # k = 11, 21, 31
                    ax.text(
                        x = 0, 
                        y = -5,
                        s= "k = ", 
                        fontsize=15, 
                        ha='left', 
                        va='bottom',
                        weight='bold'
                    )
                    # kmers read
                    ax.text(
                        x = 0, 
                        y = -10,
                        s= "kmers read", 
                        fontsize=15, 
                        ha='left', 
                        va='bottom',
                        weight='bold'
                    )

                    ax.set_xlabel("k", fontsize = 15, ha='left')
                    if 'time' in parameter:
                        ax.set_ylabel("Running Time (Seconds)", fontsize=15)
                    else:
                        ax.set_ylabel("Maximum memory constumption (kB)", fontsize=15)
                    
                    # Adjust y-axis to log scale if you notice massive differences between k values
                    # ax.set_yscale('log') 

                    # Clean up legends ?
                    # ax.legend(title="Version Run", loc="upper left")
                    # ax.get_legend().remove()

            # save
            type = ""
            if counter == 0:
                type = "sans-vs-fingerprint"
            else:
                type = "between-fingerprints"
            savefilename = f"plots/{type}_{parameter}_comparison.png"
            plt.savefig(savefilename, dpi=300)
            print("Subplots successfully generated and saved to: ", savefilename)

def plot_all_species(species):

    # Set clean visual style
    sns.set_theme(style="whitegrid")

    row = 0
    col = 0
    nrows = 3
    ncols = 1

    # comparison: sans original - fingerprint:  "FL is None or FL == 64"
    # comparison between fingerprints:  "FL is not None"
    # 'Is' nodes are not implemented
    for counter, condition in enumerate(("FL.isna() or FL == 64", "not FL.isna()")):
        for parameter in ("running time (min, adjusted)", "max memory (GB)"):

            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(
                ncols * 8, nrows * 6), sharey=False, sharex=False)

            plt.xticks((11, 21, 31))

            comparison_mode = ''
            if counter == 0:
                comparison_mode = "\n sans original vs. fingerprinting \n "
            else:
                comparison_mode = "\n between different fingerprint lengths \n "

            if "time" in parameter:
                fig.suptitle("Execution run time" + comparison_mode,
                             weight='bold', fontsize=25)
            else:
                fig.suptitle("Maximum resident set size" + comparison_mode, weight='bold', fontsize=25)\

            fig.subplots_adjust(left=0.3, right=0.9,
                                wspace=0.5, hspace=0.4, top=0.85)

            for name, counts in species.items():

                # for in-between fingerprints I now want only drosophila or ebola dataset
                if (counter == 1 and (name != "drosophila" or name != "ebola")):
                    continue

                for n_genomes in counts:
                    species_df = df.query(
                        f"({condition}) and (species == '{name}') and (n == {n_genomes})")
                    ax = axes[row]

                    print(name, "df :")
                    print(species_df)

                    if species_df.empty:
                        # ax.text(
                        #     0.5, 0.5, f"No data for {name}", ha='center', va='center')
                        # ax.set_title(name.capitalize())
                        continue

                    # Sort logically by k so the columns go left-to-right (e.g., k=11, k=21, k=31)
                    species_df = species_df.sort_values(by='k')

                    # 3. GENERATE THE SCATTER PLOT
                    # We use hue to split colors by "version" (Original vs Fingerprint)
                    sns.scatterplot(
                        data=species_df,
                        x='k',
                        y=parameter,
                        hue='version',
                        # Gives different markers (e.g., circle vs X) for extra clarity
                        style='version',
                        s=150,             # Distinct, large data points
                        ax=ax,
                        palette={'sans original': 'green',
                                 'fingerprint': 'blue',
                                 'incorrect results': 'red'},
                        edgecolor='none',
                        alpha=0.85
                    )

                    # 4. ANNOTATE DATA POINTS
                    for _, r0w in species_df.iterrows():

                        # Format annotation: display precision score or FL if applicable
                        if r0w['version'] == 'fingerprint' or r0w['version'] == 'incorrect results':
                            label_text = f"FL{int(r0w['FL'])}"
                        else:
                            label_text = ""

                        # drop text labels right next to data points
                        ax.text(
                            x=r0w['k'] + 1,
                            y=r0w[f"{parameter}"],
                            s=label_text,
                            fontsize=13,
                            ha='left',
                            va='center',
                            weight='semibold'
                        )

                    # 5. Annotate x ticks
                    tick_positions = (11, 21, 31)
                    tick_labels = []

                    # let's take the average of kmers read for each k
                    for k in (11, 21, 31):
                        column = species_df.query(f"k == {k}")['kmers read']
                        avg = column.sum() / len(column)
                        s = f"{k}\n" + format_sci(avg)
                        # the number of kmers per color set is important
                        # 64 had always correct results
                        kmers_per_cs = species_df.query(
                            # take the value from the first row
                            # square brackets!
                            f"FL == 64 and k == {k}")['kmers per color set'].iloc[0]
                        s += "\n" + format_sci(kmers_per_cs)

                        tick_labels.append(s)

                    # Set tick positions
                    ax.set_xticks(tick_positions)

                    # Make all y values be comparable from zero
                    y_max = species_df[parameter].max()
                    ax.set_ylim(bottom=-0.07*y_max, top=1.07*y_max)

                    # Set x labels and increase font size
                    ax.set_xticklabels(tick_labels, rotation=0, size=14)
                    ax.tick_params(axis='y', size=14, labelsize=16)

                    # 5. REFINING AXES AND TITLES
                    ax.set_title(f" {name}, n = {n_genomes}",
                                 fontsize=16, pad=10, weight='bold')

                    ax.set_xlabel("")

                    # x labels
                    if (col == 0):
                        ax.text(
                            x=-0.15,
                            y=-0.16,
                            s="         k\n" \
                              "kmers read\n" \
                              "avg. k-mers\n" \
                              "per color set",
                            fontsize=15,
                            ha='right',
                            va='center',
                            weight='bold',
                            transform=ax.transAxes
                        )

                        if 'time' in parameter:
                            ax.set_ylabel("Running Time \n (minutes, sum over CPUs)",
                                          fontsize=15, loc='top', labelpad=10)
                        else:
                            ax.set_ylabel(
                                "Maximum memory \n consumption (GB)", fontsize=15, loc='top')

                    # Adjust y-axis to log scale if you notice massive differences between k values
                    # ax.set_yscale('log')

                    # Clean up legends
                    if (row == 0 and col == 0):
                        ax.legend(title="Version Run", loc="best")
                    else:
                        ax.get_legend().remove()

                    # increment
                    col += 1
                    if col == ncols:
                        col = 0
                        row += 1
                    if row == nrows:
                        row = 0

            # save
            type = ""
            if counter == 0:
                type = "sans-vs-fingerprint"
            else:
                type = "between-fingerprints"
            savefilename = f"plots/{type}_{parameter}_comparison.png"
            plt.savefig(savefilename, dpi=300)
            print("Subplots successfully generated and saved to: ", savefilename)


def calculate_Lorenz_curve(raw_counts):
    # Sort values in ascending order
    arr = np.sort(raw_counts)
    n = len(arr)
    
    # Calculate the cumulative percentages
    scaled_cumulative_counts = np.cumsum(arr) / np.sum(arr)
    
    # Force the curve to anchor cleanly at the (0,0) starting origin
    lorenz_y = np.insert(scaled_cumulative_counts, 0, 0)
    lorenz_x = np.linspace(0, 1, n + 1)
    
    # Calculate the area under the Lorenz Curve using the trapezoidal rule
    area_under_lorenz = np.trapz(lorenz_y, lorenz_x)

    # Gini is the remaining proportion of the triangle area
    gini = 1.0 - 2.0 * area_under_lorenz
    
    return (lorenz_x, lorenz_y, gini)


def plot_combined_lorenz_curves(df, k_filter=None):
    """
    Plots a single combined Lorenz Curve chart where each species is represented 
    by a distinct line, allowing direct visual comparison of inequality.
    
    Parameters:
    - df: The master Pandas DataFrame containing your data.
    - k_filter: Optional integer (e.g., 11, 21, 31) if you want to isolate the 
                analysis to a specific k-mer length.
    """
    # Optional: Filter by a specific k value if desired
    if k_filter:
        plot_df = df[df['k'] == k_filter].copy()
        title_suffix = f" (k = {k_filter})"
    else:
        plot_df = df.copy()
        title_suffix = " (All Data Combined)"

    # Setup the plot environment
    plt.figure(figsize=(10, 10))
    sns.set_theme(style="whitegrid")
    
    # 1. Plot the baseline: Line of Perfect Equality (Diagonal)
    plt.plot([0, 1], [0, 1], linestyle='--', color='gray', alpha=0.7,
             label='Perfect Equality (Gini = 0.0)')
    
    # Define a distinct color palette for your species
    # (Adjust names to match your exact 'dataset' string values)
    species_list = plot_df['dataset'].unique()
    colors = sns.color_palette("Set1", n_colors=len(species_list))
    
    # 2. Iterate through each species and plot its independent curve
    for idx, species in enumerate(species_list):
        species_data = plot_df[plot_df['dataset'] == species]
        
        # Pull the raw values (Assuming column name is 'kmers per color set')
        raw_counts = species_data['kmers per color set'].dropna().values
        
        if len(raw_counts) == 0:
            continue
            
        # Sort values in ascending order
        arr = np.sort(raw_counts)
        n = len(arr)
        
        # Calculate the cumulative percentages
        scaled_cumulative_counts = np.cumsum(arr) / np.sum(arr)
        
        # Force the curve to anchor cleanly at the (0,0) starting origin
        lorenz_y = np.insert(scaled_cumulative_counts, 0, 0)
        lorenz_x = np.linspace(0, 1, n + 1)
        
        # Calculate the mathematical Gini Coefficient for this specific species
        mad = np.abs(np.subtract.outer(arr, arr)).mean()
        rmad = mad / arr.mean() if arr.mean() != 0 else 0
        gini = 0.5 * rmad
        
        # Plot the unique line for this species
        plt.plot(lorenz_x, lorenz_y, color=colors[idx], lw=2.5,
                 label=f'{species.capitalize()} (Gini = {gini:.2f})')
        
    # 3. ADVANCED TICK & GRID CUSTOMIZATION (As requested)
    # Define major ticks at every 10% (0.1) interval
    major_ticks = np.arange(0, 1.1, 0.1)
    plt.xticks(major_ticks, labels=[f"{int(i*100)}%" for i in major_ticks[0::2]], fontsize=11)
    plt.yticks(major_ticks, labels=[f"{int(i*100)}%" for i in major_ticks], fontsize=11)
    
    # Define minor ticks at every 5% (0.05) interval to catch tight curves
    minor_ticks = np.arange(0, 1.05, 0.05)
    plt.gca().set_xticks(minor_ticks, minor=True)
    plt.gca().set_yticks(minor_ticks, minor=True)
    
    # Stylize the gridlines (solid for 10% marks, subtle dotted for 5% marks)
    plt.grid(which='major', linestyle='-', linewidth=0.8, color='#e0e0e0')
    plt.grid(which='minor', linestyle=':', linewidth=0.5, color='#c0c0c0', alpha=0.7)
    
    # 4. PLOT AESTHETICS & GEOMETRY
    plt.title(f"Comparative Lorenz Curves of Species Diversity{title_suffix}", 
              fontsize=15, weight='bold', pad=15)
    plt.xlabel("Cumulative Proportion of Color Sets (Sorted by Size)", fontsize=12, labelpad=10)
    plt.ylabel("Cumulative Proportion of Total K-mers", fontsize=12, labelpad=10)
    
    plt.xlim(0, 1.0)
    plt.ylim(0, 1.0)
    
    # Place the legend cleanly inside the graph arena
    plt.legend(loc="upper left", fontsize=11, frameon=True, shadow=True)
    
    # Enforce a perfect 1:1 square aspect ratio so angles are not distorted
    plt.gca().set_aspect('equal', adjustable='box')
    
    plt.tight_layout()
    
    # Save image file
    output_filename = "combined_species_lorenz_curves.png"
    plt.savefig(output_filename, dpi=300)
    print(f"Successfully generated comparison plot! Saved as '{output_filename}'")


def add_core_kmers():
    input_folder="fingerprint_info"
    output_folder="just_columns"

    # Regex pattern: matches "Collecting core k-mers... (" followed by digits
    pattern = re.compile(r"Collecting core k-mers\.\.\.\s*\((\d+)")
    
    # Ensure input directory exists
    if not os.path.exists(input_folder):
        print(f"Error: Folder '{input_folder}' does not exist.")
        return

    # Process each matching file
    for filename in os.listdir(input_folder):
        if filename.endswith("_FL64.txt"):
            input_path = os.path.join(input_folder, filename)
            
            # Extract number from the text file
            kmer_count = None
            with open(input_path, 'r') as f:
                for line in f:
                    match = pattern.search(line)
                    if match:
                        kmer_count = match.group(1)
                        break  # Stop searching once found
            
            if kmer_count is not None:
                # Derive output filename: replace _FL64 with "" and .txt with .tsv
                tsv_filename = filename.replace("_FL64", "").replace(".txt", ".tsv")
                tsv_path = os.path.join(output_folder, tsv_filename)
                
                # Append the extracted number to the corresponding TSV file
                if os.path.exists(tsv_path):
                    with open(tsv_path, 'a') as f_out:
                        f_out.write(f"{kmer_count} 0\n")
                    print(f"Appended {kmer_count} to {tsv_path}")
                else:
                    print(f"Warning: Output file {tsv_path} does not exist.")
            else:
                print(f"Warning: Target pattern not found in {filename}")

def create_table_of_core_kmers():
    "This function collects the percentage of core kmers from each info file into a neat table"

    input_folder="fingerprint_info"

    # Regex pattern: captures the percentage digits before '%'
    # Example target line: "Collecting core k-mers... (1178562 / 60%) (15.2 min)"
    line_pattern = re.compile(r"Collecting core k-mers\.\.\..*?/\s*(\d+)%")
    
    # Regex to extract species and k value from filename 
    # Example: "ebola_n160_k31_FL64.txt" -> species="ebola", k="31"
    # Updated regex to capture: 1) species, 2) n value, 3) k value
    file_pattern = re.compile(r"^([A-Za-z0-9]+)_n(\d+)_k(\d+)_FL64\.txt$")

    records = []

    for filename in os.listdir(input_folder):
        if filename.endswith("_FL64.txt"):
            file_match = file_pattern.search(filename)
            if not file_match:
                continue

            species = file_match.group(1)
            n_val = int(file_match.group(2))  # Extracts the integer after 'n'
            k_val = int(file_match.group(3))
            input_path = os.path.join(input_folder, filename)

            # Search line by line for percentage
            with open(input_path, 'r') as f:
                for line in f:
                    match = line_pattern.search(line)
                    if match:
                        percentage = int(match.group(1))
                        records.append({
                            'species': species,
                            'n': n_val,
                            'k': k_val,
                            'percentage': percentage
                        })
                        break

    # Inside the loop:
    species = file_match.group(1)
    n_val = int(file_match.group(2))  # Extracts the integer after 'n'
    k_val = int(file_match.group(3))

    records.append({
        'species': species,
        'n': n_val,
        'k': k_val,
        'percentage': percentage
    })


    # Convert records into a raw DataFrame
    df = pd.DataFrame(records)

    # Pivot into a 2D matrix: Species on rows, k-values on columns
    # Creates a multi-index on the row axis (e.g. species on level 1, n on level 2)
    # Uses mean() by default if duplicates are found, and keeps species & n as separate columns
    pivot_df = df.pivot_table(
        index=['species', 'n'], 
        columns='k', 
        values='percentage', 
        aggfunc='first'  # Or 'mean'
    ).reset_index()

    # Clean up column axis name
    # pivot_df.columns.name = None
    
    return pivot_df
   

def extract_kmer_distribution(save_name):
    # save_name = just_columns/file_name_of_orig_results
    kmers_in_color_set = []    

    # read the columns into a list
    with open(save_name, "r") as the_two_columns:
        line = the_two_columns.readline()
        while line != "":
            l,r = map(int, line.split())
            if (l > 0):
                kmers_in_color_set.append(l)
            if (r > 0):
                kmers_in_color_set.append(r)
            line = the_two_columns.readline()
    
    return kmers_in_color_set


def plot_all_distributions(directory = "just_columns"):
    fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize=(10, 10))
    fig.suptitle("Kmer distribution into color sets (core kmers included)")
    axes = axes.flatten()
    sns.set_theme(style="whitegrid")
    fig.subplots_adjust(wspace=0.2, hspace=0.25, top=0.9)

    # For each original result plot the Lorenz curve  

    index = -1       # to which subtable we are plotting
    for s in species.keys():
        for n in species[s]:
            index += 1              # for each dataset increment
            ax = axes[index]            
            ax.set_title(f"{s}, n = {n}")

                 # Plot appearance / aesthetics settings

            major_ticks = np.arange(0, 1.1, 0.1)
            ax.set_xticks(major_ticks, labels=[f"{int(i*100)}%"
                        if round(i*10 % 2) == 0 else "" for i in major_ticks], fontsize=11)
            ax.set_yticks(major_ticks, labels=[f"{int(i*100)}%" for i in major_ticks], fontsize=11)
            
            # Define minor ticks at every 5% (0.05) interval to catch tight curves
            minor_ticks = np.arange(0, 1.05, 0.05)
            plt.gca().set_xticks(minor_ticks, minor=True)
            plt.gca().set_yticks(minor_ticks, minor=True)
            
            # Stylize the gridlines (solid for 10% marks, subtle dotted for 5% marks)
            ax.grid(which='major', linestyle='-', linewidth=0.8, color='#e0e0e0')
            ax.grid(which='minor', linestyle=':', linewidth=0.5, color='#c0c0c0', alpha=0.7)
            
            # 4. PLOT AESTHETICS & GEOMETRY
            # plt.title(f"Fix this title...", 
            #         fontsize=15, weight='bold', pad=15)
            if (index > 1):
                ax.set_xlabel("Cumulative Proportion of Color Sets (Sorted by Size)", fontsize=12, labelpad=10)
            if (index == 0 or index == 2):
                ax.set_ylabel("Cumulative Proportion of Total K-mers", fontsize=12, labelpad=10)
            
            plt.xlim(0, 1.0)
            plt.ylim(0, 1.0)
            
            # Enforce a perfect 1:1 square aspect ratio so angles are not distorted
            plt.gca().set_aspect('equal', adjustable='box')
            
            plt.tight_layout()
    

            # Plot Line of Perfect Equality (Diagonal)
            ax.plot([0, 1], [0, 1], linestyle='--', color='gray', alpha=0.7,
            label='Perfect Equality (Gini = 0.0)')
        
            for k in (11, 21, 31):

                filename =  f"{s}_n{n}_k{k}.tsv"

                # extract the kmers per color set:
                d = extract_kmer_distribution(f"{directory}/{filename}")
                lorenz_x, lorenz_y, gini = calculate_Lorenz_curve(d)

                ax.plot(lorenz_x, lorenz_y , lw=2.5,
                        label=f'k == {k}  (Gini = {gini:.2f})')

            # Create legend only after there is something plotted!
            ax.legend(loc="upper left", fontsize=11, frameon=True, shadow=True)

               

    # Save image file
    output_filename = "plots/lorenz_curves.png"
    plt.savefig(output_filename, dpi=300)
    print(f"Successfully generated comparison plot! Saved as '{output_filename}'")


def plot_all_threads(): 

    benchmark = df_prev.query(f"(version == 'sans original') and (k == 11)")
    print("printing benchmark\n:")
    print(benchmark.head())
    
    # Set clean visual style
    sns.set_theme(style="whitegrid")

    row = 0
    col = 0
    nrows = 3
    ncols = 1

    condition = "FL == 64"
    for parameter in ("running time (min, adjusted)", "max memory (GB)"):

        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(
            ncols * 8, nrows * 6), sharey=False, sharex=False)

        plt.xticks(range(33))

        comparison_mode = "\n with varying number of threads" \
        " \n for k = 11"

        if "time" in parameter:
            fig.suptitle("(Total) execution run time" + comparison_mode,
                            weight='bold', fontsize=25)
        else:
            fig.suptitle("Maximum resident set size" + comparison_mode, weight='bold', fontsize=25)\

        fig.subplots_adjust(left=0.2, right=0.9,
                            wspace=0.5, hspace=0.4, top=0.85, bottom = 0.05)

        for name, counts in species.items():

            for n_genomes in counts:
                species_df = df.query(
                    f"({condition}) and (species == '{name}') and (n == {n_genomes})")
                ax = axes[row]

                print(name, "df :")
                print(species_df)

                if species_df.empty:
                    # ax.text(
                    #     0.5, 0.5, f"No data for {name}", ha='center', va='center')
                    # ax.set_title(name.capitalize())
                    continue

                # Sort by the number of threads so the columns go left-to-right (e.g., k=11, k=21, k=31)
                species_df = species_df.sort_values(by='threads')

                # 3. GENERATE Line Plot
                # We use hue to split colors by "version" (Original vs Fingerprint)
                sns.lineplot(
                    data=species_df,
                    x='threads',
                    y=parameter,
                    hue='version',
                    # Gives different markers (e.g., circle vs X) for extra clarity
                    style='version',
                    ax=ax,
                    palette={'sans original': 'green',
                                'fingerprint': 'blue',
                                'incorrect results': 'red'},
                    linewidth=2.5,
                    marker="o",
                    markersize=8
                    )

                # Add a horizontal benchmark line obtained from the original results.
                # Musia byt ' ' !!!                 tu   a tu
                value = benchmark.loc[benchmark["species"] == name, parameter].item()   # -extract the one value
                print("VALUE =", value)
                ax.axhline(y=value, color='red', linestyle='-', linewidth=2, label='sans original, 32 threads')

                # 5. Annotate x ticks
                tick_positions = range(0, 33, 2)
                ax.set_xticks(tick_positions)

                # Make all y values be comparable from zero
                y_max = species_df[parameter].max()
                ax.set_ylim(bottom=-0.07*y_max, top=1.07*y_max)

                ax.tick_params(axis='y', size=14, labelsize=16)

                # 5. REFINING AXES AND TITLES
                ax.set_title(f" {name}, n = {n_genomes}",
                                fontsize=16, pad=10, weight='bold')

                ax.set_xlabel("concurrent threads", fontsize = 15)

                if 'time' in parameter:
                    ax.set_ylabel("Running Time \n (minutes, sum over CPUs)",
                                    fontsize=15, loc='top', labelpad=10)
                else:
                    ax.set_ylabel(
                        "Maximum memory \n consumption (GB)", fontsize=15, loc='top')

                # Clean up legends
                if (row == 0 and col == 0):
                    ax.legend(title="Version Run", loc="best")
                else:
                    ax.get_legend().remove()

                # increment
                col += 1
                if col == ncols:
                    col = 0
                    row += 1
                if row == nrows:
                    row = 0

        # save
        type = "multithreading"
        savefilename = f"plots/{type}_{parameter}.png"
        plt.savefig(savefilename, dpi=300)
        print("Subplots successfully generated and saved to: ", savefilename)


# this should only be done once!
done = True     # done indeed - ok
if not done:
    add_core_kmers()


if __name__ == "__main__":
    species = {"typhimurium":[1000, 2000], "drosophila":[12], "ebola":[160]}

    # plot_all_species(species)
    plot_all_distributions()
    # plot_all_threads()

    # Run and display the table of core kmers
    df_percentages = create_table_of_core_kmers()
    print("Percentages of core kmers table:")
    df_percentages.sort_values(by = "species")
    print(df_percentages)