import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import seaborn as sns

def ode_system(X, t, alpha, beta, gamma, delta):
    
    x, y = X
    dotx = alpha * x - beta * x * y
    doty = gamma * x * y - delta * y
    return [dotx, doty]

def simulate_ode(parms, x0, y0, t):
    alpha, beta, gamma, delta = parms
    return odeint(ode_system, [x0, y0], t, args=(alpha, beta, gamma, delta))

def main():
    # Parameters
    alpha = 0.6
    beta = 0.1
    gamma = 0.03
    delta = 0.5
    x0 = 100.0
    y0 = 40.0
    sigma = 20.0
    Time_pred = np.linspace(0, 300, 100)

    # Simulate ODE
    X = simulate_ode([alpha, beta, gamma, delta], x0, y0, Time_pred)
    df = pd.DataFrame(X, columns=['Prey', 'Predator'])
    df['Time'] = Time_pred
    # Melt the DataFrame except for the 'time' column
    df_melted = pd.melt(df, id_vars=['Time'], var_name='species', value_name='counts')

    # Add noise
    np.random.seed(123)
    df_noise = df.copy()
    df_noise['Prey'] += np.random.normal(0, sigma, len(df))
    df_noise['Predator'] += np.random.normal(0, sigma, len(df))
    df_noise_melted = pd.melt(df_noise, id_vars=['Time'], var_name='species', value_name='counts')

    
    # Plot
    sns.lineplot(data=df_melted, x='Time', y='counts', hue='species')
    sns.scatterplot(data=df_noise_melted, x='Time', y='counts', hue='species')
    plt.xlabel('Time')
    plt.ylabel('Population')
    plt.legend()
    plt.show()

    # Save data
    with open('data.pkl', 'wb') as f:
        pickle.dump(df, f)

main()


