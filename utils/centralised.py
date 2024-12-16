# -*- coding: utf-8 -*-
""" Centralised Kaplan-Meier
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter


# Read data
data_dir = os.path.join(os.getcwd(), 'vtg_km', 'local')
data_path1 = os.path.join(data_dir, 'data1.csv')
data_path2 = os.path.join(data_dir, 'data2.csv')
df1 = pd.read_csv(data_path1)
df2 = pd.read_csv(data_path2)

# Combine and filter
df = df1._append(df2, ignore_index=True)
df = df[df['COHORT_DEFINITION_ID'] == 1029]

# Create a KaplanMeierFitter object
kmf = KaplanMeierFitter()

# Fit the survival data
kmf.fit(
    list(df['TIME_AT_RISK'].values),
    event_observed=list(df['MORTALITY_FLAG'].values)
)
df_km_c = kmf.event_table

# Plot the survival curve
kmf.plot()
plt.show()
