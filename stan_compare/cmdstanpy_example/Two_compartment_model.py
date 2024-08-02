import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import seaborn as sns
from cycler import cycler

def ode_system(X, t, alpha, beta, gamma, delta):
    x, y = X
    dotx = alpha * x - beta * x
    doty = gamma * y - delta * y
    return [dotx, doty]

def simulate_ode(parms, x0, y0, t):
    alpha, beta, gamma, delta = parms
    return odeint(ode_system, [x0, y0], t, args=(alpha, beta, gamma, delta))

def main():
    # Parameters
    alpha = 0.012
    beta  = 0.018
    gamma = 0.012
    delta = 0.005
    x0 = 1500.0
    y0 = 700.0
    eps1 = 0.63
    eps2 = 0.47
    Time_pred = np.linspace(0, 300, 300)
    # sample time points from Time_pred for observed data
    Time_obs = np.random.choice(Time_pred, 50, replace=False)

    # Simulate ODE
    X = simulate_ode([alpha, beta, gamma, delta], x0, y0, Time_pred)
    df = pd.DataFrame(X, columns=['Comp1', 'Comp2'])
    df['Time'] = Time_pred

    # Add noise
    np.random.seed(7798)
    df_noise = df.copy()
    df_noise['Log_Comp1'] = np.log(df['Comp1'])
    df_noise['Log_Comp2'] = np.log(df['Comp2']) 
    df_noise['Log_Comp1'] += np.random.normal(0, eps1, len(df))
    df_noise['Log_Comp2'] += np.random.normal(0, eps2, len(df))
    df_noise['Comp1'] = np.exp(df_noise['Log_Comp1'])
    df_noise['Comp2'] = np.exp(df_noise['Log_Comp2'])

    # sample data for observed time points
    df_obs = df_noise.loc[df_noise['Time'].isin(Time_obs)]
    
    # Plot
    fig, ax = plt.subplots()

    # Define a color cycle
    colors = plt.cm.tab10.colors
    ax.set_prop_cycle(cycler('color', colors))

    for i, color in zip(range(2), colors):
        sns.scatterplot(data=df_obs, x='Time', y=f'Comp{i+1}', ax=ax, color=color)
        sns.lineplot(data=df, x='Time', y=f'Comp{i+1}', ax=ax, color=color)
        ax.plot([], [], '-o', color=color, label=f'Comp{i+1}')

    plt.title('Two-compartment model')
    plt.xlabel('Time')
    plt.ylabel('Counts')
    plt.yscale('log')
    plt.grid(alpha=0.5, linewidth=0.5)
    plt.legend()
    # Save the plot
    plt.savefig('Two_compartment_model.png')

    # Save data
    with open('data.pkl', 'wb') as f:
        pickle.dump(df, f)

main()
