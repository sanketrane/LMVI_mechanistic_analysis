import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import multiprocessing
import seaborn as sns
import time

def ode_system(X, t, alpha, beta, gamma, delta):
    x, y = X
    dotx = alpha * x - beta * x
    doty = gamma * y - delta * y
    return [dotx, doty]

def simulate_ode(parms, x0, y0, t):
    alpha, beta, gamma, delta = parms
    return odeint(ode_system, [x0, y0], t, args=(alpha, beta, gamma, delta))

def simulate_and_sample(prior_means, prior_sds, Time_pred, enter_seed=897):
    # Parameters
    alpha = np.random.normal(prior_means[0], prior_sds[0])
    beta = np.random.normal(prior_means[1], prior_sds[1])
    gamma = np.random.normal(prior_means[2], prior_sds[2])
    delta = np.random.normal(prior_means[3], prior_sds[3])
    x0 = np.random.normal(prior_means[4], prior_sds[4])
    y0 = np.random.normal(prior_means[5], prior_sds[5])
    eps1 = np.random.normal(prior_means[6], prior_sds[6])
    eps2 = np.random.normal(prior_means[7], prior_sds[7])

    # Simulate ODE
    X = simulate_ode([alpha, beta, gamma, delta], x0, y0, Time_pred)
    df = pd.DataFrame(X, columns=['Comp1', 'Comp2'])
    
    # Data transformation
    # We assume that the counts are Log-normally distributed
    df['Comp1'] = np.log(df['Comp1'])
    df['Comp2'] = np.log(df['Comp2'])

    #  Add noise -- Log-normal noise (additive in log scale)
    # set seed
    np.random.seed(enter_seed)
    df['Comp1'] += np.random.normal(0, eps1, len(df))
    df['Comp2'] += np.random.normal(0, eps2, len(df))
    
    # Rename the columns
    df['Time'] = Time_pred

    # Flatten the data
    df = pd.melt(df, id_vars='Time', var_name='Compartment', value_name='counts')

    return df

def sim(niter=10000, num_cores=4, enter_seed=897):
    # Priors -
    # We assume that the parameters follow a normal distribution with the following mean and sd
    # parameters: alpha, beta, gamma, delta, x0, y0, eps1, eps2
    prior_means = [0.012, 0.018, 0.012, 0.005, 1500.0, 700.0, 0.63, 0.47]
    prior_sds = [0.005, 0.005, 0.005, 0.005, 300.0, 100.0, 0.05, 0.05]

    Time_pred = np.linspace(0, 300, 200).astype(int)

    # Number of parallel simulations
    num_simulations = niter

    # Run simulations in parallel using multiprocessing Pool
    with multiprocessing.Pool(processes=num_cores) as pool:
        results = pool.starmap(simulate_and_sample, [(prior_means, prior_sds, Time_pred, enter_seed) for _ in range(num_simulations)])

    # Combine results into a single DataFrame
    df_parallel = pd.concat(results, ignore_index=True)
    return df_parallel


if __name__ == "__main__":
    # Run the main function and time it
    start_time = time.time()
    # Set seed
    seed = 599 # this seed made sure that 45 out of 50 time obs are unique
    np.random.seed(seed)
    
    # define number of cores for parallel processing
    num_cores = multiprocessing.cpu_count() - 2

    sims_parallel = sim(niter=10000, num_cores=num_cores, enter_seed=seed)  
    end_time = time.time()
    print(f"Execution time for main function: {end_time - start_time} seconds")

    # sample time points from sims_parallel for observed data
    Time_obs = np.random.choice(sims_parallel['Time'], 50, replace=False)
    
    # sort and check for repeated time points
    Time_obs = np.sort(Time_obs)
    Time_obs = np.unique(Time_obs)
    print(Time_obs.shape)

    # sample data for observed time points
    sims_obs = sims_parallel.loc[sims_parallel['Time'].isin(Time_obs)]

    # convert data to linear scale
    sims_obs_linear = sims_obs.copy()
    sims_obs_linear['counts'] = np.exp(sims_obs['counts'])


    # # Save the results as pickle and csv
    # sims_parallel.to_pickle('save_df/Two_compartment_model_parallel.pkl')
    # sims_parallel.to_csv('save_df/Two_compartment_model_parallel.csv', index=False)

    # Plot the results
    # Calculate mean and standard deviation
    summary_stats = sims_obs_linear.groupby(['Time', 'Compartment']).agg(
        mean_counts=('counts', 'mean'),
        std_counts=('counts', 'std')
    ).reset_index()


    sns.set_style('whitegrid')
    plt.figure(figsize=(12, 6))
    # Plot mean with error bars
    sns.scatterplot(x='Time', y='mean_counts', hue='Compartment', data=summary_stats, marker='o', s=60)
    for compartment in summary_stats['Compartment'].unique():
        subset = summary_stats[summary_stats['Compartment'] == compartment]
        plt.errorbar(subset['Time'], subset['mean_counts'], yerr=subset['std_counts'], fmt='o', label=f'{compartment} error')
    
    plt.xlabel('Time')
    plt.ylabel('Counts')
    plt.title('Two Compartment Model')
    plt.yscale('log')
    plt.legend()
    plt.savefig('plots/Two_compartment_model_parallel.png')
    # #plt.show()   

    print("Done!")