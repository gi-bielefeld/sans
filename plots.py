# Adrian Zetocha        23.6.2026

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


os.makedirs("plots", exist_ok=True)

df = pd.read_csv("TABLE.csv")
df.rename(columns={"dataset" : "species"}, inplace = True)
print(df.head())

species = {"typhimurium":[1000, 2000], "drosophila":[12], "ebola":[160]}


def format_sci(val):
    formatted = f"{val:.2e}"  # Produces e.g. '1.45e+08'
    base, exp = formatted.split('e+')
    return f"{base} x 10^{int(exp)}"


def plot_all_species(species):

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
                        ax.text(0.5, 0.5, f"No data for {species}", ha='center', va='center')
                        ax.set_title(species.capitalize())
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
    plt.xticks(major_ticks, labels=[f"{int(i*100)}%" for i in major_ticks], fontsize=11)
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
    fig, axes = plt.subplots(nrows = 1, ncols = 4, figsize=(16, 5), sharey=True, sharex=True)
    fig.suptitle("Kmer distribution into color sets")
    axes = axes.flatten()
    sns.set_theme(style="whitegrid")


    # For each original result plot the Lorenz curve  

    index = -1       # to which subtable we are plotting
    for s in species.keys():
        for n in species[s]:
            index += 1              # for each dataset increment
            ax = axes[index]            
            ax.set_title(f"{s}, n = {n}")

            
                 # Plot appearance / aesthetics settings

            major_ticks = np.arange(0, 1.1, 0.1)
            plt.xticks(major_ticks, labels=[f"{int(i*100)}%" for i in major_ticks], fontsize=11)
            plt.yticks(major_ticks, labels=[f"{int(i*100)}%" for i in major_ticks], fontsize=11)
            
            # Define minor ticks at every 5% (0.05) interval to catch tight curves
            minor_ticks = np.arange(0, 1.05, 0.05)
            plt.gca().set_xticks(minor_ticks, minor=True)
            plt.gca().set_yticks(minor_ticks, minor=True)
            
            # Stylize the gridlines (solid for 10% marks, subtle dotted for 5% marks)
            plt.grid(which='major', linestyle='-', linewidth=0.8, color='#e0e0e0')
            plt.grid(which='minor', linestyle=':', linewidth=0.5, color='#c0c0c0', alpha=0.7)
            
            # 4. PLOT AESTHETICS & GEOMETRY
            plt.title(f"Fix this title...", 
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
    

            # Plot Line of Perfect Equality (Diagonal)
            ax.plot([0, 1], [0, 1], linestyle='--', color='gray', alpha=0.7,
            label='Perfect Equality (Gini = 0.0)')
        
            for k in (11, 21, 31):

                filename =  f"{s}_n{n}_k{k}.tsv"

                # extract the kmers per color set:
                d = extract_kmer_distribution(f"{directory}/{filename}")
                lorenz_x, lorenz_y, gini = calculate_Lorenz_curve(d)

                ax.plot(lorenz_x, lorenz_y , lw=2.5,
                        label=f'{s.capitalize()} (Gini = {gini:.2f})')

               

    # Save image file
    output_filename = "plots/lorenz_curves.png"
    plt.savefig(output_filename, dpi=300)
    print(f"Successfully generated comparison plot! Saved as '{output_filename}'")

        


if __name__ == "__main__":
    species = {"typhimurium":[1000, 2000], "drosophila":[12], "ebola":[160]}
    plot_all_species(species)

    plot_all_distributions()



# Gemini help functions

# ==========================================
# 1. MOCK DATA GENERATION (To test the script)
# ==========================================
# In reality, you will build this flat dataframe by looping through your log files
# and parsing your output TSV files (to calculate precision and singletons).

# data = [
#     # Original SANS Data
#     {'Type': 'Original', 'Dataset': 'ebola', 'k': 11, 'FL': 'N/A', 'Time': 10, 'Memory_MB': 50, 'Size_MB': 5, 'Singletons_Pct': 2.1, 'Precision': 0.99, 'N': 100},
#     {'Type': 'Original', 'Dataset': 'drosophila', 'k': 11, 'FL': 'N/A', 'Time': 4500, 'Memory_MB': 12000, 'Size_MB': 3500, 'Singletons_Pct': 15.4, 'Precision': 0.98, 'N': 15},
    
#     # Fingerprint SANS Data (Ebola)
#     {'Type': 'Fingerprint', 'Dataset': 'ebola', 'k': 11, 'FL': 16, 'Time': 8, 'Memory_MB': 40, 'Size_MB': 5, 'Singletons_Pct': 2.1, 'Precision': 0.92, 'N': 100},
#     {'Type': 'Fingerprint', 'Dataset': 'ebola', 'k': 11, 'FL': 32, 'Time': 9, 'Memory_MB': 45, 'Size_MB': 5, 'Singletons_Pct': 2.1, 'Precision': 0.95, 'N': 100},
#     {'Type': 'Fingerprint', 'Dataset': 'ebola', 'k': 11, 'FL': 64, 'Time': 9.5, 'Memory_MB': 48, 'Size_MB': 5, 'Singletons_Pct': 2.1, 'Precision': 0.98, 'N': 100},
    
#     # Fingerprint SANS Data (Drosophila)
#     {'Type': 'Fingerprint', 'Dataset': 'drosophila', 'k': 11, 'FL': 16, 'Time': 3000, 'Memory_MB': 8000, 'Size_MB': 3500, 'Singletons_Pct': 15.4, 'Precision': 0.85, 'N': 15},
#     {'Type': 'Fingerprint', 'Dataset': 'drosophila', 'k': 11, 'FL': 32, 'Time': 3500, 'Memory_MB': 9000, 'Size_MB': 3500, 'Singletons_Pct': 15.4, 'Precision': 0.91, 'N': 15},
#     {'Type': 'Fingerprint', 'Dataset': 'drosophila', 'k': 11, 'FL': 64, 'Time': 4000, 'Memory_MB': 10000, 'Size_MB': 3500, 'Singletons_Pct': 15.4, 'Precision': 0.96, 'N': 15}
# ]

# df = pd.DataFrame(data)

# # Create the Index column you requested (e.g., "ebola_k11")
# df['Dataset_k'] = df['Dataset'] + '_k' + df['k'].astype(str)

# # ==========================================
# # 2. BUILDING THE REQUESTED TABLES
# # ==========================================

# print("--- Building Tables ---")

# # Filter for just Original data for Table 1a and 1b
# df_orig = df[df['Type'] == 'Original'].set_index('Dataset_k')

# # Table 1.a) Time Table (Original)
# table_1a = df_orig[['Time', 'Size_MB']].copy()
# table_1a.columns = ['Execution_Time_sec', 'Dataset_Size_MB']
# print("\nTable 1.a (Original Time & Size):\n", table_1a)

# # Table 1.b) Memory Table (Original)
# table_1b = df_orig[['Memory_MB', 'Singletons_Pct', 'Size_MB']].copy()
# print("\nTable 1.b (Original Memory & Features):\n", table_1b)


# # Filter for Fingerprint data for Table 2a and 2b
# df_fp = df[df['Type'] == 'Fingerprint']

# # Table 2.a) Time Table (Columns are FLs)
# # Using pivot_table to reshape the data so FLs become columns
# table_2a = df_fp.pivot_table(index='Dataset_k', columns='FL', values='Time')
# table_2a.columns = [f"Time_FL{col}" for col in table_2a.columns]
# print("\nTable 2.a (Fingerprint Time by FL):\n", table_2a)

# # Table 2.b) Memory Table (Columns are FLs)
# table_2b = df_fp.pivot_table(index='Dataset_k', columns='FL', values='Memory_MB')
# table_2b.columns = [f"Memory_FL{col}" for col in table_2b.columns]
# print("\nTable 2.b (Fingerprint Memory by FL):\n", table_2b)


# # ==========================================
# # 3. PLOTTING FUNCTION
# # ==========================================

# def plot_benchmark_scatter(data, x_col, y_col, hue_col='Dataset_k', 
#                            title="Benchmark Plot", x_label="", y_label="", 
#                            use_log_scale=False, filename=None):
#     """
#     Creates a scatter plot with colored dots per dataset, labels next to dots,
#     and optional log-log scaling.
#     """
#     plt.figure(figsize=(10, 7))
    
#     # Create the scatter plot using Seaborn for automatic coloring and legends
#     ax = sns.scatterplot(
#         data=data, 
#         x=x_col, 
#         y=y_col, 
#         hue=hue_col,
#         s=100,         # Size of the dots
#         palette="tab10", 
#         edgecolor="black",
#         alpha=0.8
#     )

#     # Add text labels next to each point
#     # We iterate over the dataframe to plot the dataset name near its coordinate
#     for i, row in data.iterrows():
#         plt.text(
#             x=row[x_col], 
#             y=row[y_col], 
#             s=str(row[hue_col]), 
#             fontsize=9,
#             ha='left',      # Horizontal alignment
#             va='bottom',    # Vertical alignment
#             padding=3
#         )

#     # Apply Log-Log scale if requested
#     if use_log_scale:
#         plt.xscale('log')
#         plt.yscale('log')
#         title += " (Log-Log Scale)"

#     plt.title(title, fontsize=14, pad=15)
#     plt.xlabel(x_label if x_label else x_col, fontsize=12)
#     plt.ylabel(y_label if y_label else y_col, fontsize=12)
#     plt.grid(True, which="both", ls="--", alpha=0.4)
    
#     # Move legend outside the plot so it doesn't overlap with data
#     plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
#     plt.tight_layout()

#     if filename:
#         plt.savefig(filename, dpi=300)
#         print(f"Saved plot to {filename}")
#     else:
#         plt.show()

# # ==========================================
# # 4. EXECUTING THE PLOTS
# # ==========================================

# if __name__ == "__main__":
    
#     # Plot 1: Old vs New Execution Time (Comparing against Dataset Size)
#     # We can plot the raw flat 'df' here to see Original vs Fingerprint on the same graph
#     plot_benchmark_scatter(
#         data=df,
#         x_col='Size_MB',
#         y_col='Time',
#         hue_col='Type',  # Color by Original vs Fingerprint instead of dataset
#         title="Execution Time vs Dataset Size",
#         x_label="Dataset Size (MB)",
#         y_label="Execution Time (Seconds)",
#         use_log_scale=True,  # Crucial for Ebola vs Drosophila
#         filename="time_vs_size_log.png"
#     )

#     # Plot 2: Precision vs Fingerprint Length
#     # Using only the fingerprinting data
#     plot_benchmark_scatter(
#         data=df_fp,
#         x_col='FL',
#         y_col='Precision',
#         hue_col='Dataset_k',
#         title="Precision vs Fingerprint Length",
#         x_label="Fingerprint Length (Bits)",
#         y_label="Precision Score",
#         use_log_scale=False, # Linear is better here since FL is just 16, 32, 64, 128
#         filename="precision_vs_fl_linear.png"
#     )




    

# def parse_time_string(time_str):
#     """Converts the 'm:ss' or 'h:mm:ss' output from time -v into raw seconds."""
#     parts = time_str.split(':')
#     if len(parts) == 2:
#         return int(parts[0]) * 60 + float(parts[1])
#     elif len(parts) == 3:
#         return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
#     return 0.0

# def extract_log_info(filepath):
#     """Scans a time log file and extracts memory (MB) and execution time (seconds)."""
#     mem_pattern = re.compile(r"Maximum resident set size \(kbytes\):\s+(\d+)")
#     time_pattern = re.compile(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s+([\d:.]+)")
    
#     max_memory_kb = None
#     wall_time_sec = None
    
#     with open(filepath, 'r') as f:
#         for line in f:
#             mem_match = mem_pattern.search(line)
#             if mem_match:
#                 max_memory_kb = int(mem_match.group(1))
#                 continue
                
#             time_match = time_pattern.search(line)
#             if time_match:
#                 wall_time_sec = parse_time_string(time_match.group(1))
                
#     # Convert memory to Megabytes for easier reading
#     memory_mb = max_memory_kb / 1024 if max_memory_kb else None
    
#     return memory_mb, wall_time_sec

# def build_dataframe():
#     """Iterates through both info folders and builds a flat, tidy DataFrame."""
#     data_rows = []
    
#     # 1. Parse Original SANS logs
#     for filepath in glob.glob("original_info/*.txt"):
#         filename = os.path.basename(filepath)
#         mem, time = extract_log_info(filepath)
        
#         # Example filename: typhimurium_n1000_k11.txt
#         parts = filename.replace('.txt', '').split('_')
#         dataset = parts[0]
#         n_genomes = int(parts[1].replace('n', ''))
#         k_val = int(parts[2].replace('k', ''))
        
#         data_rows.append({
#             'Type': 'Original',
#             'Dataset': dataset,
#             'N': n_genomes,
#             'k': k_val,
#             'FL': 'N/A', # Original doesn't have fingerprint length
#             'Memory_MB': mem,
#             'Time_Sec': time
#         })

#     # 2. Parse Fingerprint SANS logs
#     for filepath in glob.glob("fingerprint_info/*.txt"):
#         filename = os.path.basename(filepath)
#         mem, time = extract_log_info(filepath)
        
#         # Example filename: typhimurium_n1000_k11_FL32.txt
#         parts = filename.replace('.txt', '').split('_')
#         dataset = parts[0]
#         n_genomes = int(parts[1].replace('n', ''))
#         k_val = int(parts[2].replace('k', ''))
#         fl_val = int(parts[3].replace('FL', ''))
        
#         data_rows.append({
#             'Type': 'Fingerprint',
#             'Dataset': dataset,
#             'N': n_genomes,
#             'k': k_val,
#             'FL': str(fl_val),
#             'Memory_MB': mem,
#             'Time_Sec': time
#         })
        
#     return pd.DataFrame(data_rows)

# # ==========================================
# # EXECUTION & PLOTTING SKETCH
# # ==========================================

# if __name__ == "__main__":
    # # 1. Build the tidy dataframe
    # df = build_dataframe()
    # print("Data successfully loaded!")
    # print(df.head())

    # # 2. Plotting Memory Comparison (Original vs Fingerprint)
    # # Let's filter to just Typhimurium at k=11 to keep the graph readable
    # plot_data = df[(df['Dataset'] == 'typhimurium') & (df['k'] == 11)].copy()
    
    # # For a clean comparison, we merge the 'Type' and 'FL' columns for the legend
    # plot_data['Version'] = plot_data.apply(
    #     lambda row: 'Original' if row['Type'] == 'Original' else f"Fingerprint (FL={row['FL']})", 
    #     axis=1
    # )

    # plt.figure(figsize=(10, 6))
    
    # # Seaborn makes grouped bar charts effortless with tidy data
    # sns.barplot(
    #     data=plot_data, 
    #     x='N', 
    #     y='Memory_MB', 
    #     hue='Version'
    # )
    
    # plt.title("Memory Usage Comparison: Original vs Fingerprint SANS (Typhimurium, k=11)")
    # plt.xlabel("Number of Genomes (N)")
    # plt.ylabel("Maximum Resident Set Size (MB)")
    # plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    # plt.tight_layout()
    # plt.savefig("memory_comparison.png", dpi=300)
    # print("Plot saved as memory_comparison.png")






