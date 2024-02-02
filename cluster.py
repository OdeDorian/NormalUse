from astropy.io import fits
import numpy as np
import os
import sys
from tqdm import tqdm
import math
import Dorian
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import warnings
# from astropy.wcs import WCS
from sklearn.isotonic import IsotonicRegression
import corner
import pandas as pd
import emcee_sample

axy = Dorian.Dorian()

directory = '/Users/administrator/acorns-master/examples/output_files/output_fits'
filename = 'acorns_hierarchical_tree_38.fits'
filename_non = 'acorns_nonhierarchical_tree_110.fits'
GaiaSet = './Input/1700633435680O-result.csv'

# Load original 12CO data
hdu_12 = fits.open('./Output/Program1_12.fits')[0]
data_12 = hdu_12.data
hdr_12 = hdu_12.header

# Load original 13CO data
hdu_13 = fits.open('/Users/administrator/Python/Input/Program1_13.fits')[0]
hdr_13 = hdu_13.header
data_13 = hdu_13.data
image_file = fits.open('/Users/administrator/Python/Input/Program1_13.fits')[0]
image_hdr = image_file.header

cluster = Dorian.Dorian()

"""计算属性"""


def axis_range(hdr, V):
    axis_right = (hdr['NAXIS1'] - hdr['CRPIX1']) * hdr['CDELT1']
    axis_left = (0 - hdr['CRPIX1']) * hdr['CDELT1']
    axis_top = ((hdr['NAXIS2'] - hdr['CRPIX2']) * hdr['CDELT2'])
    axis_bottom = ((0 - hdr['CRPIX2']) * hdr['CDELT2'])
    if V:
        ext = (axis_left, axis_right, axis_top, axis_bottom)
    else:
        ext = (axis_left, axis_right, axis_bottom, axis_top)
    return ext


def v2c(header, velocity):
    velocity = np.array(velocity)
    # nc = header['NAXIS3']
    v0 = header['CRVAL3']
    c0 = header['CRPIX3']
    dv = header['CDELT3'] / 1e3
    channel = (velocity - v0) / dv + c0 - 1
    return channel


def matching(cube):
    lon, lat, v, I, dis = [], [], [], [], []
    L_center = 0
    B_center = 0
    V_center = 0
    Intensity_total = 0
    for i in range(len(cube)):
        Intensity_total += cube[i][2]
    for i in range(len(cube)):
        pixel_l = int(cube[i][0])
        pixel_b = int(cube[i][1])
        pixel_v = int(v2c(hdr_13, cube[i][4]))
        Intensity = cube[i][2] / Intensity_total
        L_center += cluster.c2p_l(hdr_13, pixel_l) * Intensity
        B_center += cluster.c2p_b(hdr_13, pixel_b) * Intensity
        V_center += cube[i][4] * Intensity
        dis.append(cube[i][4])
        lon.append(pixel_l)
        lat.append(pixel_b)
        v.append(pixel_v)
        I.append(Intensity)
    # 计算速度的加权平均值
    V_average = np.average(dis, weights=I)
    V_average = np.array(V_average)
    V_dispersion = np.sqrt(np.average((dis - V_average)**2, weights=I))
    pixel = list(zip(lat, lon))
    pixels = list(set(pixel))
    pos = list(zip(v, lat, lon))
    return pos, pixels, round(L_center, 2), round(B_center, 2), round(V_center, 2), round(V_dispersion, 2)


def CutFits(ipa, pos, ext_N):
    hdu = fits.open(ipa)
    data = hdu[0].data
    header = hdu[0].header
    data = np.squeeze(data)
    new_header = header.copy()

    CR = np.zeros(6, dtype=int)
    v, lat, lon = zip(*pos)
    CR[0] = max(1, np.min(v) - 5)
    CR[1] = min(header['NAXIS3'] - 1, np.max(v) + 5)
    CR[2] = max(1, np.min(lat) - ext_N)
    CR[3] = min(header['NAXIS2'] - 1, np.max(lat) + ext_N)
    CR[4] = max(1, np.min(lon) - ext_N)
    CR[5] = min(header['NAXIS1'] - 1, np.max(lon) + ext_N)
    sub_data = data[CR[0] - 1:CR[1], CR[2] - 1:CR[3], CR[4] - 1:CR[5]]

    new_header['NAXIS1'] = sub_data.shape[2]
    new_header['NAXIS2'] = sub_data.shape[1]
    new_header['NAXIS3'] = sub_data.shape[0]

    new_header['CRPIX1'] = header['CRPIX1'] - CR[4] + 1
    new_header['CRPIX2'] = header['CRPIX2'] - CR[2] + 1
    new_header['CRPIX3'] = header['CRPIX3'] - CR[0] + 1

    v = v - CR[0]
    lon = lon - CR[4]
    lat = lat - CR[2]
    pixel = list(zip(lat, lon))
    pos = list(zip(v, lat, lon))

    return sub_data, new_header, pixel, pos


def moment(data, hdr):
    new_hdr = hdr.copy()
    delta = new_hdr['CDELT3']
    for k in list(new_hdr.keys()):
        if k[0] == 'C' and (k[-1] == '4' or k[-1] == '3'):
            if k in new_hdr:
                del new_hdr[k]
    del new_hdr['NAXIS3']
    new_hdr['BUNIT'] = 'K km/s'
    delta_scaled = delta / 1e3
    dataOut = np.nansum(data, axis=0) * delta_scaled
    return dataOut, new_hdr


def MaxI_Map(data):
    MaxIMap = np.nanmax(data, axis=0)
    MaxIMap[np.isnan(MaxIMap)] = 0
    return MaxIMap


def excitation_temperature(data_Twelve, hdr_Twelve):
    hdr_N = hdr_Twelve.copy()
    MaxI = MaxI_Map(data_Twelve)
    TexMap = np.divide(5.532, np.log(1 + (5.532 / (MaxI + 0.819))))
    T_mean = np.nanmean(TexMap)
    T_peak = np.nanmax(TexMap)
    hdr_N['BUNIT'] = 'K'
    return TexMap, T_mean, T_peak


def optical_depth(data, TexMap, iso):
    MaxIMap = MaxI_Map(data)
    if iso:
        T = 5.28864
        fml_1 = np.divide(MaxIMap, T)
        fml_2 = np.divide(1, (np.exp(T / TexMap) - 1)) - 0.167667
        tau = - np.log(1 - np.divide(fml_1, fml_2))
    else:
        T = 5.26852
        fml_1 = np.divide(MaxIMap, T)
        fml_2 = np.divide(1, (np.exp(T / TexMap) - 1)) - 0.169119
        tau = - np.log(1 - np.divide((fml_1, fml_2)))
    Q = np.divide(tau, 1 - np.exp(-tau)) * np.divide(1, 1 - np.exp(-T / TexMap))
    return tau, Q


def column_density_13(mom, Q):
    Y = 5.0e5
    unit_ = 1e20
    cd_iso13 = 3.0e14 * Y * Q * mom
    cd_iso13unit = cd_iso13 / unit_
    cd_iso13_mean = np.nanmean(cd_iso13unit)
    cd_iso13_max = np.nanmax(cd_iso13unit)
    return cd_iso13, cd_iso13unit, cd_iso13_mean, cd_iso13_max


def column_density_12(mom):
    X = 1.8
    CDMap = np.multiply(mom, X)
    CD_max = np.nanmax(CDMap)
    CD_mean = np.nanmean(CDMap)
    # hdr['BUNIT'] = '$10^{20} cm^{-2}$'
    return CDMap, CD_mean, CD_max


def Distance(ipa_star, hdr, data, pixel, opa_c, opa_d, star, on_sig, off_sig, on_N):
    print('starting to calculate distance !')
    df = pd.read_csv(ipa_star)
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
    index_1 = np.where((Ag == np.inf) | (lon == np.inf) | (lat == np.inf) | (Ag_max == np.inf) | (Ag_min == np.inf)
                       | (d_i == np.inf) | (d_i_err == np.inf))
    for i in columns:
        columns[i] = np.delete(columns[i], index_1)
    print('Data with infinite-value has been removed, proceeding......')
    index_2 = np.where((Ag == np.nan) | (lon == np.nan) | (lat == np.nan) | (Ag_max == np.nan) | (Ag_min == np.nan)
                       | (d_i == np.nan) | (d_i_err == np.nan))
    for i in columns:
        columns[i] = np.delete(columns[i], index_2)
    print('Data with nan-value has been removed, proceeding......')
    image_hdr = hdr
    image_data = data
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
    indices_1 = np.where(columns['d'] > star)
    for i in columns:
        columns[i] = np.delete(columns[i], indices_1)
    print('Detailed selection has completed !')
    del lon, lat, Ag, Ag_max, parallax, parallax_error, d_i, d_i_err

    """从这里开始对数据更细节处理>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
    """=================================================="""
    """从轮廓上看"""
    # lat, lon = zip(*pixel)
    # l_center = int(np.mean(lon))
    # b_center = int(np.mean(lat))
    # l_ran = (image_hdr['NAXIS1'] - 1, 0)
    # b_ran = (0, image_hdr['NAXIS2'] - 1)
    #
    # sub_cube_1 = image_data[b_center:b_ran[1], l_ran[1]:l_center]
    # sub_cube_2 = image_data[b_center:b_ran[1], l_center:l_ran[0]]
    # sub_cube_3 = image_data[b_ran[0]:b_center, l_ran[1]:l_center]
    # sub_cube_4 = image_data[b_ran[0]:b_center, l_center:l_ran[0]]
    #
    # pixel = np.array(pixel)
    # pixel_1 = pixel[np.where((pixel[:, 1] < l_center) & (pixel[:, 0] >= b_center))]
    # pixel_1[:, 0] = pixel_1[:, 0] - b_center
    # pixel_2 = pixel[np.where((pixel[:, 1] >= l_center) & (pixel[:, 0] >= b_center))]
    # pixel_2[:, 0] = pixel_2[:, 0] - b_center
    # pixel_2[:, 1] = pixel_2[:, 1] - l_center
    # pixel_3 = pixel[np.where((pixel[:, 1] < l_center) & (pixel[:, 0] < b_center))]
    # pixel_4 = pixel[np.where((pixel[:, 1] >= l_center) & (pixel[:, 0] < b_center))]
    # pixel_4[:, 1] = pixel_4[:, 1] - l_center
    #
    #
    # alpha_1 = process_sub_cube(sub_cube_1, pixel_1, 150)
    # alpha_2 = process_sub_cube(sub_cube_2, pixel_2, 150)
    # alpha_3 = process_sub_cube(sub_cube_3, pixel_3, 150)
    # alpha_4 = process_sub_cube(sub_cube_4, pixel_4, 150)
    #
    # beta_1 = process_cube_alpha(alpha_1, hdr, data, 3 * rms, columns, axy, l_ran[1], b_center)
    # beta_2 = process_cube_alpha(alpha_2, hdr, data, 3 * rms, columns, axy, l_center, b_center)
    # beta_3 = process_cube_alpha(alpha_3, hdr, data, 3 * rms, columns, axy, l_ran[1], b_ran[0])
    # beta_4 = process_cube_alpha(alpha_4, hdr, data, 3 * rms, columns, axy, l_center, b_ran[0])
    #
    # P1 = process_columns(columns, beta_1[0], beta_1[1])
    # P2 = process_columns(columns, beta_2[0], beta_2[1])
    # P3 = process_columns(columns, beta_3[0], beta_3[1])
    # P4 = process_columns(columns, beta_4[0], beta_4[1])
    #
    # A_out = np.hstack([P1[0], P2[0], P3[0], P4[0]])
    # D_out = np.hstack([P1[1], P2[1], P3[1], P4[1]])
    # L_out = np.hstack([P1[2], P2[2], P3[2], P4[2]])
    # B_out = np.hstack([P1[3], P2[3], P3[3], P4[3]])
    #
    # index_on = beta_1[0] + beta_2[0] + beta_3[0] + beta_4[0]
    # print(len(index_on))
    # for i in columns:
    #     columns[i] = columns[i][index_on]
    #     columns[i] = columns[i][~np.isnan(columns[i])]
    #     print(len(columns[i]))
    """整体图"""
    lat, lon = zip(*pixel)
    alpha = np.zeros_like(data)
    alpha[lat, lon] = 1
    beta = data.copy()
    beta[lat, lon] = np.nan
    rms_off = np.sqrt(np.nanmean(beta ** 2))
    gamma = data.copy()
    index = np.where(alpha == 0)
    gamma[index] = np.nan
    rms_on = np.sqrt(np.nanmean(gamma ** 2))
    index_on, index_off = [], []
    for i in range(len(columns['b'])):
        if alpha[axy.p2c_b(hdr, columns['b'][i]), axy.p2c_l(hdr, columns['l'][i])] == 1 and \
                (data[axy.p2c_b(hdr, columns['b'][i]), axy.p2c_l(hdr, columns['l'][i])] > on_sig * rms_on):
            index_on.append(i)
        if (alpha[axy.p2c_b(hdr, columns['b'][i]), axy.p2c_l(hdr, columns['l'][i])] == 0) and \
                (data[axy.p2c_b(hdr, columns['b'][i]), axy.p2c_l(hdr, columns['l'][i])] < off_sig * rms_off):
            index_off.append(i)

    AG_delta = np.multiply(1/2, columns['Ax'][index_off] - columns['An'][index_off])
    cut_off = np.where(AG_delta >= 0.05)
    D_out = columns['d'][index_off]
    D_out = D_out[cut_off]
    A_out = columns['A'][index_off]
    A_out = A_out[cut_off]
    L_out = columns['l'][index_off]
    L_out = L_out[cut_off]
    B_out = columns['b'][index_off]
    B_out = B_out[cut_off]
    AG_delta = np.multiply(1/2, columns['Ax'][index_off] - columns['An'][index_off])[cut_off]

    ir = IsotonicRegression()
    weights = np.divide(1, AG_delta ** 2)
    x = ir.fit_transform(D_out, A_out, sample_weight=weights)

    plt.figure()
    plt.scatter(D_out, x, color='black', s=0.7)
    plt.scatter(columns['d'][index_on], columns['A'][index_on], color='red', s=0.7)

    for i in columns:
        columns[i] = columns[i][index_on]

    D = columns['d']
    A = columns['A']
    L = columns['l']
    B = columns['b']
    ax = columns['Ax']
    an = columns['An']
    der = columns['der']

    print('The number of on cloud stars is : ', len(A))
    if len(A) < on_N:
        return 0

    baseline = ir.predict(D)
    A = A - baseline
    plt.scatter(D, baseline, color='orange', s=0.7)

    zero = np.where(A > -4)
    A = A[zero]
    D = D[zero]
    L = L[zero]
    B = B[zero]
    ax = ax[zero]
    an = an[zero]
    der = der[zero]

    """计算emcee参数"""
    D_min = np.min(D)
    D_max = np.max(D)
    index_3 = np.argsort(D)[-50:]
    Ag_gt50 = A[index_3]
    mu2 = np.mean(Ag_gt50)
    std_gt50 = np.std(Ag_gt50)
    AG_delta = np.multiply(1 / 2, ax - an)

    plt.scatter(D, A, color='green', s=0.1)
    plt.xlim(0, star)
    # plt.show()

    """开始emcee过程"""

    emc = emcee_sample.emo()
    a = (D_min, D_max, 0.1, mu2, 0.5, std_gt50)     # 将mu1修改成0.1
    chain = emc.mcmc_sample(D, A, der, AG_delta, AG_delta ** 2, len(A), a, thin=10)

    """画图"""

    plt.figure(figsize=(14, 8))
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
    plt.savefig(opa_c, dpi=500)

    fig, (ax, bx) = plt.subplots(2)
    ax.set_position([0.1, 0.5, 0.8, 0.4])
    ax.imshow(image_data, extent=(l_range[1], l_range[0], b_range[0], b_range[1]), cmap='twilight',
              origin='lower')
    ax.contour(alpha, levels=[0.1], extent=(l_range[1], l_range[0], b_range[0], b_range[1])
               , colors='white', linewidths=0.8)
    ax.scatter(L, B, c='green', s=0.1)
    ax.scatter(L_out, B_out, c='blue', s=0.1)
    ax.tick_params(axis='x', direction='in')
    ax.tick_params(axis='y', direction='in')
    ax.set_title('Star distribution')
    ax.set_xlabel('Longitude (degree)')
    ax.set_ylabel('Latitude (degree)')

    bins = np.linspace(0, star, 200)
    counts, bin_edges = np.histogram(D, bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    average_y1 = []
    average_y2 = []
    for i in range(len(bin_centers)):
        y1_values = A[(D >= bin_edges[i]) & (D < bin_edges[i + 1])]
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
    bx.scatter(bin_centers, average_y1, color='green', s=0.7)
    # bx.scatter(bin_centers, average_y2, color='blue', s=0.5)
    bx.tick_params(axis='x', direction='in')
    bx.tick_params(axis='y', direction='in')
    bx.set_xlabel('Distance (pc)')
    bx.set_ylabel('$A_G (mag)$')
    x_lim = star
    x_m = np.median(chain[:, 0]).item()
    x_up = np.percentile(chain[:, 0], 84).item()
    x_lw = np.percentile(chain[:, 0], 16).item()
    mu_1 = chain[:, 1].mean()
    mu_2 = chain[:, 3].mean()
    bx.set_xlim(0, x_lim)
    bx.set_ylim(-1, 3.5)
    bx.vlines(x=x_m, ymin=-1, ymax=2.6, color='black', linestyles='solid', linewidths=0.7)
    bx.hlines(y=mu_1, xmin=0, xmax=x_m, color='blue', linestyles='dashed', linewidths=0.4)
    bx.hlines(y=mu_2, xmin=x_m, xmax=x_lim, color='blue', linestyles='dashed', linewidths=0.4)
    x_color = np.linspace(x_up, x_lw, 1000)
    bx.fill_between(x_color, -1, 2.6, color='cyan', alpha=0.8)
    bx.text(x_m, 2.6, '$%.2f^{+%.2f}_{-%.2f}$ pc' % (x_m, x_up - x_m, x_m - x_lw), ha='center', va='bottom')
    plt.savefig(opa_d, dpi=500)
    # plt.close()
    x_m = round(x_m, 1)
    print("the cloud's distance is %f pc\n----------" % x_m)
    return x_m


def mass(dist, reso, N_pixels, S_pixel, cd_iso13):
    conv_dist = 3.085678e18  # conversion between 1pc and 1cm
    distance = dist * conv_dist
    mu = 2.83
    m_H = 1.674e-27
    theta_MB = reso / 3600  # resolution in second
    area = N_pixels * S_pixel  # each pixel is 0.25 Acr-Min**2
    D = distance * math.sqrt(4 * area / 3600 - theta_MB ** 2) * (math.pi / 180)
    Radius = round((D / 2) / conv_dist, 1)
    M = mu * m_H * math.pi * (D / 2) ** 2 * cd_iso13 * 1e20
    M_sun = 1.989e30
    M_mc = M / M_sun
    M_mc = round(M_mc)
    print('Mass = %i times solar mass' % M_mc)
    return M_mc, Radius


def pipeline_properties(m, star, ext_N, on_sig, off_sig, on_N):
    sys.stdout = open(os.devnull, 'w')
    Gaia_Set = './Input/1700633435680O-result.csv'
    print('Making Table  >>>>>>>>>>>')
    warnings.filterwarnings("ignore")
    df = pd.DataFrame(columns=['ID', 'Area (SquareMin)', 'L_center', 'B_center', 'Vcenter (km/s)', 'Radius(pc)',
                               'Distance (pc)', 'Mass_12 (Solar Mass)', 'Mass_13 (Solar Mass)',
                               'N_mean_H2 (1e20/SquareCM)', 'N_max_H2 (1e20/SquareCM)', 'Tex_Mean (K)',
                               'Tex_Peak (K)', 'Velocity_dispersion (km/s)'])
    IDSet = [np.nan] * m
    DistanceSet = [np.nan] * m
    MassSet = [np.nan] * m
    TpeakSet = [np.nan] * m
    CDmeanSet = [np.nan] * m
    CDmaxSet = [np.nan] * m
    AreaSet = [np.nan] * m
    TmeanSet = [np.nan] * m
    RadiusSet = [np.nan] * m
    LcenterSet = [np.nan] * m
    BcenterSet = [np.nan] * m
    VcenterSet = [np.nan] * m
    VdispersionSet = [np.nan] * m
    Mass_12Set = [np.nan] * m

    directory = '/Users/administrator/acorns-master/examples/output_files/output_fits'
    df_book = pd.read_csv('~/Documents/book1.csv')
    id = df_book.loc[:, 'Index'].values
    a = 0
    region_1 = [16, 21, 93, 38, 85, 49, 9]
    region_2 = [107, 149, 61, 98, 75, 42, 57, 193, 94]      # 1 rms
    region_3 = [52, 5, 6, 33, 2, 32, 53, 213, 13, 11, 81, 48]       # 1 rms
    region_total = region_1 + region_2 + region_3
    for k in tqdm(region_total):
        filename = 'acorns_hierarchical_tree_%i.fits' % k
        cube = fits.open(os.path.join(directory, filename))[1].data
        match = matching(cube)
        cut13 = CutFits('/Users/administrator/Python/Output/Program1_13.fits', pos=match[0], ext_N=ext_N)
        mom = moment(cut13[0], cut13[1])
        """这里的剪切的12数据的速度范围有问题"""
        cut12 = CutFits('/Users/administrator/Python/Output/Program1_12.fits', pos=match[0], ext_N=ext_N)
        mom_12 = moment(cut12[0], cut12[1])
        cd_12 = column_density_12(mom_12[0])
        mass_12 = mass(189, 50, len(cut13[2]), 0.25, cd_12[1])
        tex = excitation_temperature(cut12[0], cut12[1])
        tau = optical_depth(cut13[0], tex[0], iso=True)
        cd = column_density_13(mom[0], tau[1])
        # dist = Distance(Gaia_Set, mom[1], mom[0], cut13[2], './Output/Fi/%iCM.pdf' % k,
        #                 './Output/Fi/%iD.pdf' % k, star, on_sig, off_sig, on_N)
        Mass = mass(204, 50, len(cut13[2]), 0.25, cd[2])

        IDSet[a] = k
        DistanceSet[a] = int(204)
        MassSet[a] = Mass[0]
        LcenterSet[a] = match[2]
        BcenterSet[a] = match[3]
        VcenterSet[a] = match[4]
        TpeakSet[a] = np.round(tex[2], 2)
        TmeanSet[a] = np.round(tex[1], 2)
        CDmeanSet[a] = np.round(cd[2], 2)
        CDmaxSet[a] = np.round(cd[3], 2)
        RadiusSet[a] = Mass[1]
        AreaSet[a] = len(match[1]) * 0.25
        VdispersionSet[a] = match[5]
        Mass_12Set[a] = mass_12[0]

        data = pd.Series \
            ({
                'ID': IDSet[a],
                'Area (SquareMin)': AreaSet[a],
                'L_center': LcenterSet[a],
                'B_center': BcenterSet[a],
                'Vcenter (km/s)': VcenterSet[a],
                'Radius(pc)': RadiusSet[a],
                'Distance (pc)': DistanceSet[a],
                'Mass_13 (Solar Mass)': MassSet[a],
                'Mass_12 (Solar Mass)': Mass_12Set[a],
                'N_mean_H2 (1e20/SquareCM)': CDmeanSet[a],
                'N_max_H2 (1e20/SquareCM)': CDmaxSet[a],
                'Tex_Mean (K)': TmeanSet[a],
                'Tex_Peak (K)': TpeakSet[a],
                'Velocity_dispersion (km/s)': VdispersionSet[a]
            })
        df = df._append(data, ignore_index=True)
        df.to_csv('./MC_Catalogue_Settled.csv', index=True)
        a += 1
    else:
        pass
    sys.stdout = sys.__stdout__
    print('Property Table has been established.')


def pipeline_plot(catalogue):
    directory = '/Users/administrator/acorns-master/examples/output_files/output_fits'
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 6)
    ext = axis_range(hdr_13, V=False)
    index = 0
    df = pd.read_csv(catalogue)
    lis = np.round(df.loc[:, 'ID'].values)
    longitude = df.loc[:, 'L_center'].values
    latitude = df.loc[:, 'B_center'].values
    dist = df.loc[:, 'Distance (pc)'].values
    mom = moment(data_13, hdr_13)
    ClrMap = ('dimgrey', 'cyan', 'pink', 'lime', 'yellow', 'orange')
    ext_cm = ClrMap * 7
    # latus = np.flipud(mom[0])
    # img = ax.imshow(latus, extent=(ext[0], ext[1], ext[2], ext[3]), cmap='Purples', vmin=0, vmax=7)
    latus = np.empty_like(mom[0])
    for i in tqdm(range(len(lis))):
        file_name = 'acorns_hierarchical_tree_%i.fits' % lis[i]
        cube = fits.open(os.path.join(directory, file_name))[1].data
        match = matching(cube)
        lat, lon = zip(*match[1])
        patus = np.zeros_like(mom[0])
        patus[lat, lon] = 1
        latus[lat, lon] += mom[0][lat, lon]

        """画轮廓线 + 输出代号"""
        if (dist[i] < 300) and (dist[i] > 0):
            ax.contour(patus, levels=1, extent=(ext[0], ext[1], ext[2], ext[3]), colors=ext_cm[i], linewidths=0.4,
                       alpha=0.5)
            ax.annotate('%i' % lis[i], xy=(longitude[i], latitude[i]), xytext=(longitude[i], latitude[i]), fontsize=8)
        elif dist[i] < 600:
            ax.contour(patus, levels=[0.01], extent=(ext[0], ext[1], ext[2], ext[3]), colors='orange', linewidths=0.5)
            # ax.annotate('%i' % lis[i], xy=(longitude[i], latitude[i]), xytext=(longitude[i], latitude[i]))
        elif dist[i] < 900:
            ax.contour(patus, levels=[0.01], extent=(ext[0], ext[1], ext[2], ext[3]), colors='green', linewidths=0.5)
            # ax.annotate('%i' % lis[i], xy=(longitude[i], latitude[i]), xytext=(longitude[i], latitude[i]))
        elif dist[i] < 1200:
            ax.contour(patus, levels=[0.01], extent=(ext[0], ext[1], ext[2], ext[3]), colors='cyan', linewidths=0.5)
            # ax.annotate('%i' % lis[i], xy=(longitude[i], latitude[i]), xytext=(longitude[i], latitude[i]))
        index += 1
    # ax.invert_xaxis()
    latus[np.where(latus == 0)] = np.nan
    rms = np.sqrt(np.nanmean(latus**2))
    latus = np.flipud(latus)
    img = ax.imshow(latus, extent=(ext[0], ext[1], ext[2], ext[3]), cmap='Purples', vmin=0, vmax=5.5)
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%d°'))
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%d°'))
    ax.set_xlabel("Galactic Longitude")
    ax.set_ylabel("Galactic Latitude")
    divider = make_axes_locatable(ax)
    cax1 = divider.append_axes("right", size="5%", pad=0)
    cax1 = fig.colorbar(img, cax=cax1)
    cax1.set_label('K km/s')
    fig.savefig('./MC_Map.png', dpi=500)
    fig.savefig('./MC_Map.pdf')
    fig.show()


def MaskMatch(sub_cube, hdr, pos):
    v, lat, lon = zip(*pos)
    # mask = np.where(sub_cube != pos, np.nan, data_o)
    mask = np.empty_like(sub_cube)
    mask[v, lat, lon] = sub_cube[v, lat, lon]
    mask[(np.where(mask == 0))] = np.nan
    hdu_n = fits.PrimaryHDU(mask, header=hdr)
    hdu_n.writeto('./Region_3_13.fits', overwrite=True)
    print('match complete !')

# pipeline_properties(100, 1400, 100, 1, 1.5, 0)
pipeline_plot('./MC_Catalogue_Settled.csv')


# region_1 = [16, 21, 93, 38, 85, 49, 9]      # 1 rms
# region_2 = [107, 149, 61, 98, 75, 42, 57, 193, 94]      # 1 rms
# region_3 = [52, 5, 6, 33, 2, 32, 53, 213, 13, 11, 81, 48]       # 1 rms
# region_4 = [69, 23, 3]      # 0 rms
# region_5 = [26, 17]     # 1 rms
# region_6 = [15, 12]     # 0 rms
# region_7 = [30]
# region_8 = [95]
# region_total = [16, 21, 93, 38, 85, 49, 9, 107, 149, 61, 98, 75, 42, 57, 193, 94, 52,
#                 5, 6, 33, 2, 32, 53, 213, 13, 11, 81, 48]
# directory = '/Users/administrator/acorns-master/examples/output_files/output_fits'
# v, lat, lon = [], [], []
# for i in region_3:
#     filename = 'acorns_hierarchical_tree_%i.fits' % i
#     cube = fits.open(os.path.join(directory, filename))[1].data
#     match = matching(cube)
#     v_, lon_, lat_ = zip(*match[0])
#     v += v_
#     lon += lon_
#     lat += lat_
# pos_total = list(zip(v, lon, lat))
# hdu = CutFits('/Users/administrator/Python/Output/Program1_13.fits', pos_total, 100)
# MaskMatch(hdu[0], hdu[1], hdu[3])
# mom = moment(hdu[0], hdu[1])
# dist = Distance(GaiaSet, mom[1], mom[0], hdu[2], './region_8_C.pdf', './region_8.pdf', 1400, 1, 1, 0)

"""选择更加细节的off星区域"""


def process_sub_cube(sub_cube, pixel, num):
    cube_alpha = np.zeros_like(sub_cube)
    for i in range(len(sub_cube[0])):
        gamma = pixel[pixel[:, 0] == i, 1]
        if len(gamma) > 0:
            max_pl = np.max(gamma)
            min_pl = np.min(gamma)
            max_out = min(max_pl + num, len(cube_alpha[1]) - 1)
            min_out = max(min_pl - num, 0)
            if not max_pl + 5 > max_out or min_pl - 5 < min_out:
                cube_alpha[i, min_out:min_pl - 5] = 1
                cube_alpha[i, max_pl + 5:max_out] = 1
            elif max_pl + 5 < max_out:
                cube_alpha[i, max_pl + 5:max_out] = 1
            elif min_pl - 5 > min_out:
                cube_alpha[i, min_out:min_pl - 5] = 1
            cube_alpha[i, min_pl:max_pl] = 2
    return cube_alpha


def process_cube_alpha(cube_alpha, hdr, data, Nrms, columns, axy, l_offset, b_offset):
    """找到on和off星的星表索引"""
    hdr_N = hdr.copy()
    hdr_N['CRPIX1'] = hdr['CRPIX1'] - l_offset + 1
    hdr_N['CRPIX2'] = hdr['CRPIX2'] - b_offset + 1
    hdr_N['NAXIS1'] = cube_alpha.shape[1] + 1
    hdr_N['NAXIS2'] = cube_alpha.shape[0] + 1
    index_on, index_off = [], []
    for i in range(len(columns['b'])):
        if columns['l'][i] < axy.c2p_l(hdr_N, 0) and columns['l'][i] > axy.c2p_l(hdr_N, hdr_N['NAXIS1'] - 2) and columns['b'][i] > axy.c2p_b(hdr_N, 0) and columns['b'][i] < axy.c2p_b(hdr_N, hdr_N['NAXIS2'] - 2):
            if cube_alpha[axy.p2c_b(hdr_N, columns['b'][i]), axy.p2c_l(hdr_N, columns['l'][i])] == 2:
                index_on.append(i)
            if (cube_alpha[axy.p2c_b(hdr_N, columns['b'][i]), axy.p2c_l(hdr_N, columns['l'][i])] == 1) and (data[axy.p2c_b(hdr, columns['b'][i]), axy.p2c_l(hdr, columns['l'][i])] < Nrms):
                index_off.append(i)
    return index_on, index_off


def process_columns(columns, index_on, index_off):
    D_out = columns['d'][index_off]
    A_out = columns['A'][index_off]
    Lon_out = columns['l'][index_off]
    Lat_out = columns['b'][index_off]
    AG_delta = np.multiply(1/2, columns['Ax'][index_off] - columns['An'][index_off])

    ir = IsotonicRegression()
    weights = np.divide(1, AG_delta ** 2)
    ir.fit_transform(D_out, A_out, sample_weight=weights)

    baseline = ir.predict(columns['d'][index_on])
    columns['A'][index_on] = columns['A'][index_on] - baseline
    return A_out, D_out, Lon_out, Lat_out


"""修正距离后的质量"""


# mass(189, 50, 1763*4, 0.25, 3.3)
