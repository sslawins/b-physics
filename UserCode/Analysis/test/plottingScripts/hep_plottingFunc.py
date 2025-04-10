import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import uproot as upr
import awkward as ak
from matplotlib.colors import LogNorm
import shutil
import os
import mplhep as hep
import matplotlib.pyplot as plt
import numba as nb
import re



hep.style.use("CMS")
params = {'legend.fontsize': 'x-large',
              'figure.figsize': (10, 7),
              'axes.labelsize': 'x-large',
              'axes.titlesize':'x-large',
              'xtick.labelsize':'x-large',
              'ytick.labelsize':'x-large'}
plt.rcParams.update(params)


def sanitize_filename(filename):
    return re.sub(r'\W+', '_', filename)

def histogram_1D(data, column, bins, xlabel, ylabel, title, fig_path, save=False):
    plt.figure(figsize=(20, 15))
    plt.hist(data[column], bins=bins, histtype='step', color='blue')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
    hep.cms.text("Private", fontsize=40)

    if save:
        sanitized_title = sanitize_filename(title)
        plt.savefig(os.path.join(fig_path, sanitized_title + '.png'))
    
    plt.show()

    
def histogram_2D(data, column1, column2, bins, xlabel, ylabel, title, fig_path, save=False, log_scale=False):
    plt.figure(figsize=(20, 15))
    if log_scale:
        h = plt.hist2d(data[column1], data[column2], bins=bins, norm=LogNorm())
        plt.colorbar(h[3], ax=plt.gca())
    else:
        plt.hist2d(data[column1], data[column2], bins=bins)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
    hep.cms.text("Private", fontsize=40)

    if save:
        sanitized_title = sanitize_filename(title)
        plt.savefig(os.path.join(fig_path, sanitized_title + '.png'))
    
    plt.show()

def histo1D_prompt_displaced(data_prompt, data_displaced, column, bins, xlabel, ylabel, title, fig_path, save=False):
    plt.figure(figsize=(20, 15))
    plt.hist(data_prompt[column], bins=bins, histtype='step', color='blue', label='Prompt')
    plt.hist(data_displaced[column], bins=bins, histtype='step', color='red', label='Displaced')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    
    hep.cms.text("Private", fontsize=40)

    if save:
        sanitized_title = sanitize_filename(title)
        plt.savefig(os.path.join(fig_path, sanitized_title + '.png'))
    
    plt.show()