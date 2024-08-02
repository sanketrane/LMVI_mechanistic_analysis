import pandas as pd
import numpy as np
import os
from cmdstanpy import CmdStanModel
import pickle

# time sequence for predictions 
ts_pred = np.linspace(4, 30, num=500)
numPred = len(ts_pred)

# Read the CSV file
data_df = pd.read_csv("artf_obs_noisy.csv")

# Create a dictionary containing the data
data_file = {
  'numObs1': len(data_df['Times']),
  'solve_time': data_df['Times'],
  'CAR_MZ_counts': data_df['MZ'],
  'CAR_GC_counts': data_df['GC'],
  'CAR_MZN2_counts':  data_df['MZ-N2'],
  'CAR_GCN2_counts':  data_df['GC'],
  'ts_pred': ts_pred,
  'numPred': numPred
}

# write the stan model code to a file