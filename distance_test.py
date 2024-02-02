import sys
import time
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import math
from spectral_cube import SpectralCube
import astropy.units as u
import os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker
import matplotlib.patches as patches
from astropy.table import Table
import corner
import pandas as pd
import emcee_sample
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from tqdm import tqdm
import warnings
import Dorian
from scipy.ndimage import zoom
from distfit import distfit
from fitter import Fitter
import numpy as np
import matplotlib.pyplot as plt
from sklearn.isotonic import IsotonicRegression


axy = Dorian.Dorian()

def process_sub_cube(sub_cube, l_offset, b_offset, image_hdr, columns):
    alpha, beta = np.where(sub_cube == 1)
    theta = np.array(list(zip(alpha, beta)))
    alpha = np.unique(alpha)

    index_alpha = []
    for k in alpha:
        gamma = theta[theta[:, 0] == k, 1] + l_offset

        max_pl = np.max(gamma) + 100
        min_pl = np.min(gamma) - 100

        max_l = axy.c2p_l(image_hdr, min_pl)
        min_l = axy.c2p_l(image_hdr, max_pl)

        max_b = axy.c2p_b(image_hdr, k + b_offset) + image_hdr['CDELT2']
        min_b = axy.c2p_b(image_hdr, k + b_offset)

        index_alpha += list(np.where(
            (columns['b'] > min_b) & (columns['b'] < max_b) & (columns['l'] < max_l) & (columns['l'] > min_l)))


    index_beta = []
    for sublist in index_alpha:
        index_beta.extend(sublist)
    index_beta = list(set(index_beta))

    return index_beta


def process_columns(columns, index_beta, image_data, image_hdr, axy):
    ag = columns['A'][index_beta]
    delta_ag = np.multiply(1/2, columns['Ax'][index_beta] - columns['An'][index_beta])
    dist = columns['d'][index_beta]
    lon = columns['l'][index_beta]
    lat = columns['b'][index_beta]
    der = columns['der']
    index = np.argsort(dist)
    arr = (np.array(lon)[index], np.array(lat)[index], np.array(dist)[index], np.array(ag)[index], np.array(delta_ag)[index], np.array(der)[index])
    mask = np.zeros_like(arr[0], dtype=int)
    mat = np.zeros_like(arr[0], dtype=int)
    for i in range(len(arr[0])):
        if not image_data[axy.p2c_b(image_hdr, arr[1][i]), axy.p2c_l(image_hdr, arr[0][i])] > 0:
            mask[i] = 1
        else:
            mat[i] = 1

    D_out = [d for i, d in enumerate(arr[2]) if mask[i]]
    A_out = [a for i, a in enumerate(arr[3]) if mask[i]]
    A_delta = [a for i, a in enumerate(arr[4]) if mask[i]]
    L_out = [d for i, d in enumerate(arr[0]) if mask[i]]
    B_out = [a for i, a in enumerate(arr[1]) if mask[i]]

    A_delta = np.array(A_delta)
    cuto_ff = np.where(A_delta <= 0.05)
    D_out = np.delete(D_out, cuto_ff)
    A_out = np.delete(A_out, cuto_ff)
    L_out = np.delete(L_out, cuto_ff)
    B_out = np.delete(B_out, cuto_ff)
    A_delta = np.delete(A_delta, cuto_ff)
    A_delta = np.array(A_delta)

    D_on = [d for i, d in enumerate(arr[2]) if mat[i]]
    A_on = [a for i, a in enumerate(arr[3]) if mat[i]]
    delta_ag_on = [k for i, k in enumerate(arr[4]) if mat[i]]
    L_on = [a for i, a in enumerate(arr[0]) if mat[i]]
    B_on = [a for i, a in enumerate(arr[1]) if mat[i]]
    der_on = [a for i, a in enumerate(arr[5]) if mat[i]]

    ir = IsotonicRegression()
    weights = np.divide(1, A_delta**2)
    ir.fit_transform(D_out, A_out, sample_weight=weights)

    baseline = ir.predict(D_on)
    A_on = np.array(A_on)
    D_on = np.array(D_on)
    delta_ag_on = np.array(delta_ag_on)
    der_on = np.array(der_on)
    A_on = A_on - baseline
    return A_on, D_on, delta_ag_on, L_on, B_on, A_out, D_out, L_out, B_out, der_on


def Distance(inpathS, inpathM, outpathc, outpathd):

    """
    Instruction
        inpathS: the path of gaia star catalog (csv)
        inpathM: the path of input moment fits (fits)
        outpathof: the path of the distribution of binned on & off cloud stars
        outpaths: the path of stars in moment figure
        outpathc: the path of corner map
        outpathd: the path of distance figure
    return: distance of the sample molecular cloud
    """

    """筛选数据"""

    print('starting to calculate distance !')
    df = pd.read_csv(inpathS)
    lon = df.loc[:, 'l'].values
    lat = df.loc[:, 'b'].values
    Ag = df.loc[:, 'ag_gspphot'].values
    Ag_max = df.loc[:, 'ag_gspphot_upper'].values
    Ag_min = df.loc[:, 'ag_gspphot_lower'].values
    parallax = df.loc[:, 'parallax'].values
    d_i = 1 / parallax * 1e3
    parallax_error = df.loc[:, 'parallax_error'].values
    d_i_err = 1 / parallax_error
    columns = {'l': lon, 'b': lat, 'A': Ag, 'Ax': Ag_max, 'An': Ag_min, 'd': d_i, 'der': d_i_err}
    for i in columns:
        if not isinstance(columns[i], np.ndarray):
            print('Error the input data (%s) is not a numpy ndarray' % i)
            sys.exit()
    print('Raw data columns are all arrays, proceeding......')
    index_0 = np.where((Ag == 0) | (lon == 0) | (lat == 0) | (Ag_max == 0) | (Ag_min == 0) | (parallax == 0) | (
            parallax_error == 0))
    for i in columns:
        columns[i] = np.delete(columns[i], index_0)
    print('Data with 0-value has been removed, proceeding......')
    index_1 = np.where(
        (Ag == np.inf) | (lon == np.inf) | (lat == np.inf) | (Ag_max == np.inf) | (Ag_min == np.inf) | (
                d_i == np.inf) | (d_i_err == np.inf))
    for i in columns:
        columns[i] = np.delete(columns[i], index_1)
    print('Data with infinite-value has been removed, proceeding......')
    index_2 = np.where(
        (Ag == np.nan) | (lon == np.nan) | (lat == np.nan) | (Ag_max == np.nan) | (Ag_min == np.nan) | (
                d_i == np.nan) | (d_i_err == np.nan))
    for i in columns:
        columns[i] = np.delete(columns[i], index_2)
    print('Data with nan-value has been removed, proceeding......')
    image_file = fits.open(inpathM)[0]
    image_hdr = image_file.header
    image_data = image_file.data
    l_range = (axy.c2p_l(image_hdr, image_hdr['NAXIS1'] - 1), axy.c2p_l(image_hdr, 0))
    b_range = (axy.c2p_b(image_hdr, 0), axy.c2p_b(image_hdr, image_hdr['NAXIS2'] - 1))
    print('l range =', l_range)
    print('b range =', b_range)
    indice_1 = np.where((columns['l'] < l_range[0]) | (columns['l'] > l_range[1]))
    for i in columns:
        columns[i] = np.delete(columns[i], indice_1)
    indice_2 = np.where((columns['b'] < b_range[0]) | (columns['b'] > b_range[1]))
    for i in columns:
        columns[i] = np.delete(columns[i], indice_2)
    indices_1 = np.where(columns['d'] > 2000)
    for i in columns:
        columns[i] = np.delete(columns[i], indices_1)
    print('Detailed selection has completed !')
    del lon, lat, Ag, Ag_max, parallax, parallax_error, d_i, d_i_err

    lantern = np.where(~np.isnan(image_data))
    lran = (np.min(lantern[1]), np.max(lantern[1]))
    lcenter = int(np.mean(lantern[1]))
    bran = (np.min(np.min(lantern[0]) - 10, 0), np.max(np.max(lantern[0]) + 10), image_hdr['NAXIS2'] - 1)
    bcenter = int(np.mean(lantern[0]))

    frame = np.zeros_like(image_data)
    frame[lantern] = 1
    sub_cube_1 = frame[bcenter:bran[1], lran[0]:lcenter]
    sub_cube_2 = frame[bcenter:bran[1], lcenter:lran[1]]
    sub_cube_3 = frame[bran[0]:bcenter, lran[0]:lcenter]
    sub_cube_4 = frame[bran[0]:bcenter, lcenter:lran[1]]

    index_beta_1 = process_sub_cube(sub_cube_1, lran[0], bcenter, image_hdr, columns)
    index_beta_2 = process_sub_cube(sub_cube_2, lcenter, bcenter, image_hdr, columns)
    index_beta_3 = process_sub_cube(sub_cube_3, lran[0], bran[0], image_hdr, columns)
    index_beta_4 = process_sub_cube(sub_cube_4, lcenter, bran[0], image_hdr, columns)

    P1 = process_columns(columns, index_beta_1, image_data, image_hdr, axy)
    P2 = process_columns(columns, index_beta_2, image_data, image_hdr, axy)
    P3 = process_columns(columns, index_beta_3, image_data, image_hdr, axy)
    P4 = process_columns(columns, index_beta_4, image_data, image_hdr, axy)

    A_on = np.hstack([P1[0], P2[0], P3[0], P4[0]])
    D_on = np.hstack([P1[1], P2[1], P3[1], P4[1]])
    delt_Ag = np.hstack([P1[2], P2[2], P3[2], P4[2]])
    L_on = np.hstack([P1[3], P2[3], P3[3], P4[3]])
    B_on = np.hstack([P1[4], P2[4], P3[4], P4[4]])
    A_out = np.hstack([P1[5], P2[5], P3[5], P4[5]])
    D_out = np.hstack([P1[6], P2[6], P3[6], P4[6]])
    L_out = np.hstack([P1[7], P2[7], P3[7], P4[7]])
    B_out = np.hstack([P1[8], P2[8], P3[8], P4[8]])
    der_on = np.hstack([P1[9], P2[9], P3[9], P4[9]])

    non = np.where(~np.isnan(A_on))
    A_on = A_on[non]
    D_on = D_on[non]
    L_on = L_on[non]
    B_on = B_on[non]
    der_on = der_on[non]
    delt_Ag = delt_Ag[non]

    print('The number of on cloud stars is : ', len(A_on))

    """计算emcee参数"""
    D_min = np.min(D_on)
    D_max = np.max(D_on)
    index_3 = np.argsort(D_on)[-50:]
    Ag_gt50 = A_on[index_3]
    mu2 = np.mean(Ag_gt50)
    std_gt50 = np.std(Ag_gt50)

    """开始emcee过程"""

    emc = emcee_sample.emo()
    a = (D_min + 200, D_max, 0.5, mu2, 0.5, std_gt50)
    chain = emc.mcmc_sample(D_on, A_on, der_on, delt_Ag, delt_Ag ** 2, len(A_on), a, thin=10)

    """画图"""

    plt.figure(figsize=[14, 8])
    figcor = corner.corner(chain, show_titles=True, title_kwargs={"fontsize": 10})

    axcor = np.array(figcor.axes).reshape([5, 5])

    axtil0 = 'D = ' + axcor[0, 0].title.get_text() + ' pc'
    axtil1 = '$\mu_1$ = ' + axcor[1, 1].title.get_text() + ' mag'
    axtil2 = '$\sigma_1$ = ' + axcor[2, 2].title.get_text() + ' mag'
    axtil3 = '$\mu_2$ = ' + axcor[3, 3].title.get_text() + ' mag'
    axtil4 = '$\sigma_2$ = ' + axcor[4, 4].title.get_text() + ' mag'

    axcor[0, 0].set_title(axtil0, fontsize=10)
    axcor[1, 1].set_title(axtil1, fontsize=10)
    axcor[2, 2].set_title(axtil2, fontsize=10)
    axcor[3, 3].set_title(axtil3, fontsize=10)
    axcor[4, 4].set_title(axtil4, fontsize=10)

    axcor[0, 0].axvline(chain[:, 0].mean(), linestyle='--', lw=1.0, color='Black')
    axcor[1, 1].axvline(chain[:, 1].mean(), linestyle='--', lw=1.0, color='Black')
    axcor[2, 2].axvline(chain[:, 2].mean(), linestyle='--', lw=1.0, color='Black')
    axcor[3, 3].axvline(chain[:, 3].mean(), linestyle='--', lw=1.0, color='Black')
    axcor[4, 4].axvline(chain[:, 4].mean(), linestyle='--', lw=1.0, color='Black')
    plt.savefig(outpathc, dpi=500)

    fig, (ax, bx) = plt.subplots(2)
    ax.set_position([0.1, 0.5, 0.8, 0.4])
    ax.imshow(image_data, extent=[l_range[1], l_range[0], b_range[0], b_range[1]], cmap='twilight',
              origin='lower')
    ax.scatter(L_on, B_on, color='green', s=0.1)
    ax.scatter(L_out, B_out, color='blue', s=0.1)
    ax.tick_params(axis='x', direction='in')
    ax.tick_params(axis='y', direction='in')
    ax.set_title('Star distribution')
    ax.set_xlabel('Longitude (degree)')
    ax.set_ylabel('Latitude (degree)')

    bins = np.linspace(0, 2500, 100)
    counts, bin_edges = np.histogram(columns['d'], bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    average_y1 = []
    average_y2 = []
    for i in range(len(bin_centers)):
        y1_values = A_on[(D_on >= bin_edges[i]) & (D_on < bin_edges[i + 1])]
        y2_values = A_out[(D_out >= bin_edges[i]) & (D_out < bin_edges[i + 1])]
        if y1_values.size > 0:
            average_y1.append(np.mean(y1_values))
        else:
            average_y1.append(np.nan)
        if y2_values.size > 0:
            average_y2.append(np.mean(y2_values))
        else:
            average_y2.append(np.nan)
    bx.set_position([0.1, 0.1, 0.8, 0.3])
    bx.scatter(bin_centers, average_y1, color='green', s=0.5)
    bx.scatter(bin_centers, average_y2, color='blue', s=0.5)

    bx.scatter(D_out, A_out, color='pink', s=0.5)
    bx.scatter(D_on, D_out, color='orange', s=0.5)

    bx.tick_params(axis='x', direction='in')
    bx.tick_params(axis='y', direction='in')
    bx.set_xlabel('Distance (pc)')
    bx.set_ylabel('$A_G (mag)$')
    x_lim = 2500
    x_m = np.median(chain[:, 0])
    x_up = np.percentile(chain[:, 0], 84)
    x_lw = np.percentile(chain[:, 0], 16)
    mu_1 = chain[:, 1].mean()
    mu_2 = chain[:, 3].mean()
    bx.set_xlim(0, x_lim)
    bx.set_ylim(-1, 3.5)
    bx.vlines(x=x_m, ymin=0, ymax=3, color='black', linestyles='solid')
    bx.hlines(y=mu_1, xmin=0, xmax=x_m, color='blue', linestyles='dashed')
    bx.hlines(y=mu_2, xmin=x_m, xmax=x_lim, color='blue', linestyles='dashed')
    x_color = np.linspace(x_up, x_lw, 1000)
    bx.fill_between(x_color, 0, 3, color='yellow', alpha=0.5)
    bx.text(x_m, 3, '$%.2f^{+%.2f}_{-%.2f}$ pc' % (x_m, x_up - x_m, x_m - x_lw), ha='center', va='bottom')
    plt.savefig(outpathd, dpi=500)
    x_m = round(x_m, 1)
    print("the cloud's distance is %f pc\n----------" % x_m)
    return x_m



GaiaSet = './Input/1700633435680O-result.csv'
# axy.Distance(GaiaSet, './13+match_moment.fits', './%iasdf.pdf' % 1, './%ifff.pdf' % 1)
Distance(GaiaSet, './13+match_moment.fits', './%iasdf.pdf' % 1, './%ifff.pdf' % 1)
