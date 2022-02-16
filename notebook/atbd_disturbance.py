'''
atbd_disturbance.py

NISAR Algorithm Theoretical Baseline Document (ATBD) 
FOREST DISTURBANCE ALGORITHM

Author: Josef Kellndorfer, PhD
Date: 04-Feb-2022

Description

This python module contains code elements for the ATBD notebook NISAR_Disturbance_ATBD.ipynb

This module should be imported with:


import sys
libpath=<path to directory containing this module>
sys.path.append(libpath)
from atbd_disturbance import *
'''
# Turn of warnings
# import logging
# logging.getLogger("param.main").setLevel(logging.CRITICAL)
import warnings
warnings.filterwarnings("ignore")

# Standard modules 
import os,sys
import xarray as xr
import numpy as np
import pandas as pd
import datetime as dt
import hvplot.xarray
import panel as pn
import holoviews as hv
import fsspec

from osgeo import gdal

import importlib


# User defined functions
def time_label_from_idx(ds,idx):
    label=''
    for i in idx:
        t = ds.isel({'time':i})
        l =f'{t.time.values}'.split('T')[0]
        label += f'{l} '
    return label

'''
equation template for markdown cell

\begin{align}
\label{eq1}\tag{1}
\end{align}

'''

def single_point_CP_vis(Xpoint,x,y,polarization='hv'):
    X = Xpoint.sel(x=x,y=y,method='nearest')
    Xmean_maxmin = ( X.max(dim='time') + X.min(dim='time') ) / 2.
    R = X - Xmean_maxmin
    Rm = X- X.mean()

    S = R.cumsum(dim='time')
    Sm= Rm.cumsum(dim='time')
    S.name='Cumulative Sum'
    R.name='Residuals'
    Sm.name='Cumulative Sum'
    Rm.name='Residuals'
    X.name=f'Backscatter {polarization}'

    CP = S.time.isel(time = S.argmax('time').values).values
    
    hvp = (
     X.hvplot(title='X',frame_width=250,label='gamma0') * 
     hv.HLine(float(Xmean_maxmin.values),label='min(X)+(max(X)-min(X))/2').opts(color='red') * 
     hv.HLine(float(X.mean().values),label='Mean').opts(color='gray',line_dash='dashed') + 
        R.hvplot(title='Residuals',frame_width=250)*Rm.hvplot(title='Residuals',frame_width=250).opts(color='gray',line_dash='dashed')*hv.HLine(0).opts(color='red')*hv.VLine(CP).opts(color='black') + 
     S.hvplot(title='S Curve (Cumulative Sums)',frame_width=250)*Sm.hvplot(title='S Curve (Cumulative Sums)',frame_width=250).opts(color='gray',line_dash='dashed')*hv.VLine(CP).opts(color='black')
    ).opts(title=f'Time Series of ALOS-2 Observations at Lon/lat {x}/{y}')
    
    return hvp
        