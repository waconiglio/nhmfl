import pandas as pd
import numpy as np

def integrate_bdot(bdot, area, dt):
    # integrate
    field = np.cumsum(np.asarray(bdot)/area*dt)
    # correct ending field to zero
    field -= field[-1] * np.asarray(range(len(field))) / len(field)
    return field

def read_hld_pulse(filename, columns=None, dt=1e-6, t0=-0.005, integrate_column='Bdot', area=None):
    df = pd.read_csv(filename, sep="\t", header=None, skip_blank_lines=True, comment='#', names=columns, dtype=np.float64)
    df['Time'] = df.index*dt + t0
    df.set_index('Time', inplace=True)

    # integrate field if so indicated
    if integrate_column is not None and integrate_column != '' and area is not None:
        df['Field'] = integrate_bdot(df[integrate_column].values, area, dt)
        
    return df

