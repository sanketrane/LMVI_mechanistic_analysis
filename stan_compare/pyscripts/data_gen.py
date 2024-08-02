import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Read the CSV file
df = pd.read_csv("artf_obs_data.csv")

# Set the seed
np.random.seed(4563)

# Generate CAR counts with added noise
CAR_GC_counts = np.exp(np.log(df['X0']) + np.random.normal(0, 0.58, df.shape[0]))
CAR_MZ_counts = np.exp(np.log(df['X1']) + np.random.normal(0, 0.6, df.shape[0]))
CAR_MZN2_counts = np.exp(np.log(df['X2']) + np.random.normal(0, 0.49, df.shape[0]))
solve_times = np.arange(4, 30, 0.5)

# Create a new DataFrame
new_df = pd.DataFrame({
    "Times": solve_times,
    "GC": CAR_GC_counts,
    "MZ": CAR_MZ_counts,
    "MZ-N2": CAR_MZN2_counts
})

# Melt the DataFrame for plotting
melted_df = new_df.melt(id_vars=["Times"], var_name="pop", value_name="counts")

my_palette = {
    'GC': 'tab:blue',
    'MZ': 'tab:green',
    'MZ-N2': 'tab:orange',
}

# # Plotting
# plt.figure(figsize=(10, 6))
# sns.lineplot(data=df, x=4 + df.index / 2, y='X0', color='tab:blue', legend=False)
# sns.lineplot(data=df, x=4 + df.index / 2, y='X1', color='tab:green', legend=False)
# sns.lineplot(data=df, x=4 + df.index / 2, y='X2', color='tab:orange', legend=False)
# scatter = sns.scatterplot(data=melted_df, x="Times", y="counts", hue="pop", palette=my_palette)
# plt.yscale('log')
# # Get the legend object and set the title to an empty string
# legend = scatter.legend_
# legend.set_title('')
# # Position the legend
# plt.legend(loc='upper left')
# # Set axes labels
# plt.xlabel('Time (Hours)')
# plt.ylabel('Counts')
# plt.show()
# # Write the new DataFrame to a CSV file
# new_df.to_csv("artf_obs_noisy.csv", index=False)
