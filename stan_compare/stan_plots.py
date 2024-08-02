import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from itertools import islice

# Load the stanfit object
with open('fit_file.pkl', 'rb') as f:
    fit = pickle.load(f)

# Get the list of all Stan variables
stan_variables = fit.stan_variables()

# List of all variable names in stan_variables
all_vars = list(stan_variables.keys())

# Find the index of 'sigma3' in the list
par_end_index = all_vars.index('sigma3') + 1 if 'sigma3' in all_vars else len(all_vars)

# Select variables up to 'sigma3'
selected_vars = all_vars[:par_end_index]

# Create a new dictionary with only the selected variables
selected_params = {var: stan_variables[var] for var in selected_vars}

# Calculate the number of rows and columns for the grid
n = len(selected_params)
cols = 4
rows = np.ceil(n / cols).astype(int)

# Create a figure and a grid of subplots
fig, axs = plt.subplots(rows, cols, figsize=(15, rows*3))

# Flatten the array of axes, so we can easily iterate over it
axs = axs.flatten()

# Iterate over each variable and plot the histogram
for i, var in enumerate(selected_params):
    values = selected_params[var]
    axs[i].hist(values, bins=30, edgecolor='black')
    axs[i].set_title(f'Histogram of {var}')
    axs[i].set_xlabel(var)
    axs[i].set_ylabel('Frequency')

# Remove empty subplots
if n % cols != 0:
    for ax in axs[n:]:
        fig.delaxes(ax)

plt.tight_layout()
#plt.show()

# Save the plot
saveparpath = 'plots/param_hists.png'
plt.savefig(saveparpath)

## plot predictions

# time sequence for predictions 
ts_pred = np.linspace(4, 30, num=500)
numPred = len(ts_pred)

# Read observed data from the CSV file
noisy_df = pd.read_csv("artf_obs_noisy.csv").melt(id_vars=["Times"], var_name="popln", value_name="counts")

my_palette = {
    'GC': 'tab:blue',
    'MZ': 'tab:orange',
    'MZ-N2': 'tab:green',
}

# List of unique populations
populations = noisy_df['popln'].unique()


# Rename the variables for better readability
stan_variables = {
    'GC': stan_variables['y1_mean_pred'],
    'MZ': stan_variables['y2_mean_pred'],
    'MZ-N2': stan_variables['y4_mean_pred']
    }

# List of data variables
data_vars = ['GC', 'MZ', 'MZ-N2']

# Create a grid of subplots
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(18, 4.5))
#axs = axs.flatten()

for var in data_vars:
    # Extract the relevant data for each variable
    if var in stan_variables:
        var_data = stan_variables[var]
    else:
        raise KeyError(f"The key '{var}' is not found in the fit object.")
    
    # Create a DataFrame from the extracted data
    var_pred = pd.DataFrame(var_data).melt(var_name="key", value_name="value").groupby("key").agg(
        mean_pred=pd.NamedAgg(column="value", aggfunc="mean"),
        lower_pred=pd.NamedAgg(column="value", aggfunc=lambda x: np.percentile(x, 2.5)),
        upper_pred=pd.NamedAgg(column="value", aggfunc=lambda x: np.percentile(x, 97.5))
    )

    # Plot the observed data in 3 separate subplots
    ax = axs[data_vars.index(var)]
    
    # Filter the noisy_df based on the current population
    filtered_noisy_df = noisy_df[noisy_df['popln'] == var]
    
    # Plot the noisy data
    sns.scatterplot(ax=ax, data=filtered_noisy_df, x="Times", y="counts", hue="popln", palette=my_palette)
    
    # Plot the predicted data
    axs[data_vars.index(var)].plot(ts_pred, var_pred['mean_pred'], color = my_palette[var] , label='Predictions')
    axs[data_vars.index(var)].fill_between(ts_pred, var_pred['lower_pred'], var_pred['upper_pred'], color = my_palette[var], alpha=0.3)
    #axs[data_vars.index(var)].scatter(noisy_df['Times'], noisy_df['counts'], color='red', label='Observed')
    axs[data_vars.index(var)].set_title(f'Rep+ {var} cells')
    axs[data_vars.index(var)].set_xlabel('Time since immunization (days)')
    axs[data_vars.index(var)].set_ylabel('Counts')
    axs[data_vars.index(var)].legend()
    axs[data_vars.index(var)].set_yscale('log')

    # dont show legend for noisy data
    handles, labels = axs[data_vars.index(var)].get_legend_handles_labels()
    axs[data_vars.index(var)].legend(handles=handles[1:], labels=labels[1:])



# Save the plot
savepredpath = 'plots/predictions.png'
plt.savefig(savepredpath)

# Adjust the layout
plt.tight_layout()
plt.show()


print(f"Plots saved to {saveparpath} and {savepredpath}")
print("Donzel!")
