import numpy as np
import matplotlib.pyplot as plt

# the Kernel function
def K(x):
    return np.exp(-x**2/2)/np.sqrt(2*np.pi)

# dummy dataset
dataset = np.array([1.33, 0.3, 0.97, 1.1, 0.1, 1.4, 0.4])

# x-value range for plotting KDEs
x_range = np.linspace(dataset.min()-0.3, dataset.max()+0.3, num=600)

# bandwith values for experimentation
H = [0.3, 0.1, 0.03]
n_samples = dataset.size

# line properties for different bandwith values
color_list = ['goldenrod', 'black', 'maroon']
alpha_list = [0.8, 1, 0.8]
width_list = [1.7,2.5,1.7]

plt.figure(figsize=(10,4))
# iterate over bandwith values
for h, color, alpha, width in zip(H, color_list, alpha_list, width_list):
    total_sum = 0
    # iterate over datapoints
    for i, xi in enumerate(dataset):
        total_sum += K((x_range - xi) / h)
        plt.annotate(r'$x_{}$'.format(i+1),
                     xy=[xi, 0.13],
                     horizontalalignment='center',
                     fontsize=18,
                    )
    y_range = total_sum/(h*n_samples)
    plt.plot(x_range, y_range, 
             color=color, alpha=alpha, linewidth=width, 
             label=f'{h}')

    plt.plot(dataset, np.zeros_like(dataset) , 's', 
             markersize=8, color='black')

plt.xlabel('$x$', fontsize=22)
plt.ylabel('$f(x)$', fontsize=22, rotation='horizontal', labelpad=20)
plt.legend(fontsize=14, shadow=True, title='$h$', title_fontsize=16)
plt.show()

# ---------------------------
# 
# """"""
# ---------------------------


# import seaborn as sns

# def plot_distance_distributions(atoms):
#     """
#     Compute histograms and density curves for each nucleotide type based on pairwise distances.

#     Parameters
#     ----------
#     atoms : list
#         List of atom entries, where each entry contains chain, residue ID, and 3D coordinates.

#     Returns
#     -------
#     None
#     """

#     # Compute all distances
#     distances = residue_distances(atoms)

#     # Prepare a dictionary for each nucleotide type
#     nucleotide_distances = {nuc: [] for nuc in ("A", "U", "G", "C")}

#     for res_i, res_j, d in distances:
#         # Assign distance to both nucleotides (res_i and res_j)
#         if res_i in nucleotide_distances:
#             nucleotide_distances[res_i].append(d)
#         if res_j in nucleotide_distances:
#             nucleotide_distances[res_j].append(d)

#     # Plot histograms + KDEs
#     plt.figure(figsize=(10, 6))
#     colors = {"A": "red", "U": "blue", "G": "green", "C": "orange"}

#     for nuc, dist_list in nucleotide_distances.items():
#         if len(dist_list) == 0:
#             continue

#         sns.histplot(dist_list, bins=num_bins, kde=False, stat="density",
#                      color=colors[nuc], alpha=0.3, label=f"{nuc} histogram")

#         sns.kdeplot(dist_list, color=colors[nuc], lw=2, label=f"{nuc} KDE")

#     plt.title("Distance distributions per nucleotide")
#     plt.xlabel("Distance (Å)")
#     plt.ylabel("Density")
#     plt.legend()
#     plt.show()
