import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import Dorian
from astropy.io import fits
import numpy as np
import pandas as pd
import os
import sys
import Dorian
from tqdm import tqdm
import warnings


Sundries = Dorian.Dorian()
# Sundries.AverageSpectra('./Output/Program1_13.fits', './averagemap.png')
# Sundries.MomentPlot('./my_slice.fits', 1, './pvextractor.png', contour=False)
# ran = (300, 427, 652, 1618, 1, 1471)
# Sundries.CutFits(*ran, Dimension=3, inpath='/Users/administrator/PMO/Icarus/DataBase/G140_170/DataPartI/12PartI.fits', outpath='./Output/Program1_12.fits', Layer=False)
# ran = (1, 2, 3, 4, 5, 6)
# Sundries.CutFits(*ran, Dimension=3, inpath='/Users/administrator/PMO/Dorian/314802.fits', outpath='./314802_WithoutNan_rms.fits', Layer=True)
# ran = (294, 1478, 1018, 2782)
# Sundries.CutFits(*ran, Dimension=2, inpath='/Users/administrator/PMO/Icarus/DataBase/G140_170/OverallData/12CO/G140_170_U_rms.fits', outpath='./314802_WithoutNan_rms.fits', Layer=False)
# Sundries.StackBump('/Users/administrator/PMO/Work/Work_Scripts/DBScan/G140_170_UdbscanS2P4Con1_Clean.fits', '/Users/administrator/PMO/Icarus/DataBase/G140_170/OverallData/13CO/G140_170_L.fits', '/Users/administrator/PMO/Icarus/DataBase/G140_170/OverallData/13CO/G140_170_L_rms.fits', 314802, './a.csv')
# Sundries.MomentPlot('./Overview_bvmap.fits', 2, './overview_bvmap.png', contour=False)
# ran = (278, 433, 294, 1478, 1018, 2782)
# Sundries.CutFits(*ran, Dimension=3, inpath='/Users/administrator/PMO/Icarus/DataBase/G140_170/OverallData/12CO/G140_170_U.fits', outpath='./suit_withoutNan.fits', Layer=False)
# Sundries.ChannelMap('./314802_suit.fits', './', 10, 5, autoLayout=True)
# ran = (1, 4, 2, 4, 2, 4)
# Sundries.CutFits(*ran, Dimension=3, inpath='./314802.fits', outpath='./314802_suit.fits', Layer=True)
# Sundries.MomentPlot('./314802_lvmap.fits', 1, './314802_lvmap.png', contour=False)
# hdr = fits.open('/Users/administrator/Daedalus/MC/MWISPData/ThirteenData/Daedalus13.fits')[0].header
# a = Sundries.v2c(hdr, 2)
# print(a)
# Sundries.ChannelMap('/Users/administrator/Daedalus/MC/MWISPData/TwelveData/Daedalus12.fits', '/Users/administrator/Daedalus/MC/MWISPData/TwelveData/', 3, 1.5, autoLayout=False)
# ran = (0, 20)
# mean = Sundries.CalculateTex('/Users/administrator/Daedalus/MC/MWISPData/TwelveData/Daedalus12.fits', './Tex_12.fits', FITS=True)[0]
# peak = Sundries.CalculateTex('/Users/administrator/Daedalus/MC/MWISPData/TwelveData/Daedalus12.fits', './Tex_12.fits', FITS=True)[1]
# Sundries.CalculateMoment(outfits + 'Icarus12.fits', 0, ran, outfits + 'IcarIIII.fits', ty=False)
# Sundries.MomentPlot('./Tex_12.fits', 0, './Tex.png', contour=True)
# mean = np.zeros(11)
# max = np.zeros(11)
# for i in range(3, 6):
#     Sundries.MomentPlot('./%i' % i + 'ColumnDensity.fits', 0, './%i' % i + 'CD.png', contour=True)
    # mean[i] = Sundries.Caculate12CD('./%ith_moment.fits' % i, './%i' % i, Fits=True)[0]
    # max[i] = Sundries.Caculate12CD('./%ith_moment.fits' % i, './%i' % i, Fits=True)[1]
# hdr = fits.open(inpath + '12PartI.fits')[0].header
# cutrange = (363, 406, Sundries.p2c_b(hdr, 1.5), Sundries.p2c_b(hdr, 6.5), Sundries.p2c_l(hdr, 157), Sundries.p2c_l(hdr, 147))
# cutrange = (Sundries.p2c_b(hdr, 1.5), Sundries.p2c_b(hdr, 6.5), Sundries.p2c_l(hdr, 157), Sundries.p2c_l(hdr, 147))
# Sundries.CutFits(*cutrange, Dimension=3, inpath=inpath + '18PartI.fits', outpath='/Users/administrator/Daedalus/MC/MWISPData/ThirteenData/Daedalus18.fits', Layer=False)
# Sundries.CutFits(*cutrange, Dimension=2, inpath=inpath + '18PartI_rms.fits', outpath='/Users/administrator/Daedalus/MC/MWISPData/ThirteenData/Daedalus18_rms.fits', Layer=False)
# df = pd.read_csv('/Users/administrator/PMO/Icarus/DataReduction/DataInput/Layers.csv')
# idx = df.loc[:, '_idx']
# a = []
# for i in tqdm(range(len(idx))):
#     indices = Sundries.StackBump(
#         '/Users/administrator/PMO/Icarus/DataReduction/DataInput/12PartIdbscanS2P4Con1_Clean.fits',
#         '/Users/administrator/PMO/Icarus/DataBase/G140_170/DataPartI/18PartI.fits',
#         '/Users/administrator/PMO/Icarus/DataBase/G140_170/DataPartI/18PartI_rms.fits', idx[i], outtables)
#     if indices is not None:
#         a.append(indices)
# df = pd.DataFrame(a)
# # 将DataFrame输出为CSV文件
# df.to_csv('./output18.csv', index=False)

# Sundries.ReadFit('/Users/administrator/PMO/Work/Work_Scripts/DBScan/314802_WithoutNandbscanS3P4Con1_Clean.fit')
# Sundries.MaskMatch('/Users/administrator/Daedalus/MC/MWISPData/314802_WithoutNan.fits', '/Users/administrator/PMO/Work/Work_Scripts/DBScan/314802_WithoutNandbscanS3P4Con1_Clean.fits', 1809, './1809.fits')
# Sundries.CalculateMoment('./1809.fits', 0, (-10, 20), './1809_moment.fits', ty=False)
# Sundries.MomentPlot('/Users/administrator/PMO/Dorian/1809_moment.fits', 0, './1809.png', contour=False)
# Arange12 = (1, 505, 1, 1618, 1252, 3661)
# Brange12 = (1, 505, 1, 1269, 1, 1251)
# Crange12 = (1, 1618, 1252, 3661)
# Arange13 = (1, 484, 1, 1618, 1252, 3661)
# Brange13 = (1, 484, 1, 1269, 1, 1251)
# Crange13 = (1, 1618, 1252, 3661)
# Arange18 = (1, 482, 1, 1618, 1252, 3661)
# Brange18 = (1, 482, 1, 1269, 1, 1251)

# Sundries.CalculateMoment(inpath + '18PartI.fits', 0, (-60, 20), outpath=outfits + '18PartI.fits', ty=False)

# Sundries.CalculateMoment('/Users/administrator/Daedalus/MC/MWISPData/ThirteenData/Daedalus13.fits', 0, (0, 10), '/Users/administrator/Daedalus/MC/MWISPData/ThirteenData/Moment13_0.fits', ty=False)
# Sundries.MomentPlot('/Users/administrator/Daedalus/MC/MWISPData/ThirteenData/Moment13_0.fits', 0, '/Users/administrator/Daedalus/Background/13COMoment_0.pdf', contour=False)
# Crange = (127, 674, 470, 1117)
# rmscut = (1, 1618, 1252, 3661)
# cut = (1, 482, 1, 1618, 1252, 3661)
# Sundries.CutFits(*cut, Dimension=3, inpath='/Users/administrator/PMO/Icarus/DataBase/G140_170/OverallData/C18O/G140_170_L2.fits', outpath='/Users/administrator/PMO/Icarus/DataBase/G140_170/DataPartI/18PartI.fits', Layer=False)
# Sundries.Distance(inpath + '1699597630438O-result.csv', outfits + '35998_Moment.fits', outfigures + 'star.png', outfigures + 'starsonmap.png', outfigures + 'corner.png', outfigures + 'distance.png')
# hdu = fits.open('/Users/administrator/PMO/Icarus/DataBase/G140_170/OverallData/12CO/G140_170_U.fits')[0]
# hdr = hdu.header
# Cr = (1, 505, Sundries.p2c_b(hdr, -5), Sundries.p2c_b(hdr, 5), Sundries.p2c_l(hdr, 170), Sundries.p2c_l(hdr, 140))
# Sundries.CutFits(*Cr, Dimension=3, inpath='/Users/administrator/PMO/Icarus/DataBase/G140_170/OverallData/12CO/G140_170_U.fits', outpath=outfits + 'HIIMap.fits', Layer=False)
# Sundries.Mass(700, 300, 5)





# current_script_path = os.path.realpath(__file__)
# parent_directory = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
# outfits = parent_directory + '/Icarus/DataReduction/DataOutput/Fits/'
# outfigures = parent_directory + '/Icarus/DataReduction/DataOutput/Figures/'
# outtable = parent_directory + '/Icarus/DataReduction/DataOutput/Tables/'
# inpath = parent_directory + '/Icarus/DataReduction/DataInput/'
# GaiaSet = '1700633435680O-result.csv'
#
# df = pd.read_csv(outtable + 'CloudCatalog.csv')
# L = df.loc[:, 'Distance (pc)'].values
# A = df.loc[:, 'Area (SquareMin)'].values
# N = df.loc[:, 'NMean (1e20/SquareCM)'].values
# N = N * 2 / 1.8
# R = np.empty_like(L)
# Mass = np.empty_like(L)
# for i in tqdm(range(0, len(L))):
#     Mass[i] = Icarus.Mass(L[i], A[i], N[i])[0]
#     R[i] = Icarus.Mass(L[i], A[i], N[i])[1]

# 定义文件夹路径
# import os
# import glob
# import numpy as np
# from tqdm import tqdm

# # 定义文件夹路径
# folder_path = '/Users/administrator/PMO/Icarus/DataBase/Clustered/Fits/Moment'
#
# # 使用glob模块获取所有文件路径，按文件名排序
# file_paths = sorted(glob.glob(os.path.join(folder_path, '*')))

# # 初始化mean和max为NumPy数组
# mean = np.array([])
# max = np.array([])

# # 按顺序读取每个文件
# for file_path in tqdm(file_paths):
#     result = Icarus.Caculate12CD(file_path, Fits=False)
#     mean = np.append(mean, result[0])
#     max = np.append(max, result[1])
# Sundries.ReadFit('/Users/administrator/PMO/Work/Work_Scripts/DBScan/G140_170_UdbscanS2P4Con1_Clean.fit')



# 画DBscan聚类云的轮廓图的方法
# -----------------------------------------
# Sundries.ReadFit('/Users/administrator/PMO/Work/Work_Scripts/DBScan/314802_WithoutNandbscanS3P5Con1_Clean.fit')
# image_file = fits.open('/Users/administrator/PMO/Work/Work_Scripts/DBScan/314802_WithoutNandbscanS3P5Con1_Clean.fits')[0]
# image_hdr = image_file.header
# data = image_file.data
# p_start = (image_hdr['NAXIS1'] - image_hdr['CRPIX1']) * image_hdr['CDELT1']
# p_end = (0 - image_hdr['CRPIX1']) * image_hdr['CDELT1']
# v_start = ((image_hdr['NAXIS2'] - image_hdr['CRPIX2']) * image_hdr['CDELT2'])
# v_end = ((0 - image_hdr['CRPIX2']) * image_hdr['CDELT2'])
# p_r = (p_start, p_end)
# v_r = (v_start, v_end)
#
# df = pd.read_csv('/Users/administrator/PMO/Work/Work_Scripts/DBScan/DBscanCatalog.csv')
# id = df.loc[:, '_idx'].values
# area = df.loc[:, 'area_exact'].values
# standard = np.where(area > 50)
# id = id[standard]
#
# fig, ax = plt.subplots()
# fig.set_size_inches(10, 5)
# clr = ['red', 'orange', 'yellow', 'green', 'cyan', 'blue', 'purple']
# extended_clr = clr * (1000 // len(clr)) + clr[:1000 % len(clr)]
# x = np.arange(p_r[1], p_r[0])
# y = np.arange(v_r[1], v_r[0])
# np.meshgrid(x, y)
#
# for i in tqdm(range(len(id))):
#     patus = np.zeros((image_hdr['NAXIS2'], image_hdr['NAXIS1']), dtype=int)
#     indices = np.where(data == id[i])
#     P = list(zip(indices[1], indices[2]))
#     index = list(set(tuple(i) for i in P))
#     for j in index:
#         patus[j[0]][j[1]] = 1
#     ax.contour(patus, levels=[0.1], extent=[p_r[1], p_r[0], v_r[1], v_r[0]], colors=extended_clr[i], alpha=0.6, linewidths=0.5)
#
# ax.invert_xaxis()
# ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%d°'))
# ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%d°'))
# ax.set_xlabel("Galactic Longitude")
# ax.set_ylabel("Galactic Latitude")
# os.remove('/Users/administrator/PMO/Work/Work_Scripts/DBScan/314802_WithoutNandbscanS3P5Con1_Clean.fit')
# os.remove('/Users/administrator/PMO/Work/Work_Scripts/DBScan/314802_WithoutNandbscanS3P5Con1_Clean.fits')
# os.remove('/Users/administrator/PMO/Work/Work_Scripts/DBScan/314802_WithoutNandbscanS3P5Con1.fits')
# os.remove('/Users/administrator/PMO/Work/Work_Scripts/DBScan/DBscanCatalog.csv')
# fig.savefig('./CCM.png', dpi=500)
# fig.savefig('./CCM.pdf')
# fig.show()
# -----------------------------------------

# 拼接红外数据
# -----------------------------------------
# r00data = fits.open('/Users/administrator/Daedalus/SF/WISE22/WISE22_147_6.fits')[0].data
# r01data = fits.open('/Users/administrator/Daedalus/SF/WISE22/WISE22_151_6.fits')[0].data
# r02data = fits.open('/Users/administrator/Daedalus/SF/WISE22/WISE22_155_6.fits')[0].data
# r03data = fits.open('/Users/administrator/Daedalus/SF/WISE22/WISE22_159_6.fits')[0].data
# r04data = fits.open('/Users/administrator/Daedalus/SF/WISE22/WISE22_163_6.fits')[0].data
# hdr = fits.open('/Users/administrator/Daedalus/SF/WISE22/WISE22_147_6.fits')[0].header
# rtotal01 = np.hstack((r01data, r00data))
# rtotal02 = np.hstack((r02data, rtotal01))
# rtotal03 = np.hstack((r03data, rtotal02))
# rtotal04 = np.hstack((r04data, rtotal03))
#
# r10data = fits.open('/Users/administrator/Daedalus/SF/WISE22/WISE22_147_2.fits')[0].data
# r11data = fits.open('/Users/administrator/Daedalus/SF/WISE22/WISE22_151_2.fits')[0].data
# r12data = fits.open('/Users/administrator/Daedalus/SF/WISE22/WISE22_155_2.fits')[0].data
# r13data = fits.open('/Users/administrator/Daedalus/SF/WISE22/WISE22_159_2.fits')[0].data
# r14data = fits.open('/Users/administrator/Daedalus/SF/WISE22/WISE22_163_2.fits')[0].data
# rtotal11 = np.hstack((r11data, r10data))
# rtotal12 = np.hstack((r12data, rtotal11))
# rtotal13 = np.hstack((r13data, rtotal12))
# rtotal14 = np.hstack((r14data, rtotal13))
#
# r20data = fits.open('/Users/administrator/Daedalus/SF/WISE22/WISE22_147_-2.fits')[0].data
# r21data = fits.open('/Users/administrator/Daedalus/SF/WISE22/WISE22_151_-2.fits')[0].data
# r22data = fits.open('/Users/administrator/Daedalus/SF/WISE22/WISE22_155_-2.fits')[0].data
# r23data = fits.open('/Users/administrator/Daedalus/SF/WISE22/WISE22_159_-2.fits')[0].data
# r24data = fits.open('/Users/administrator/Daedalus/SF/WISE22/WISE22_163_-2.fits')[0].data
# rtotal21 = np.hstack((r21data, r20data))
# rtotal22 = np.hstack((r22data, rtotal21))
# rtotal23 = np.hstack((r23data, rtotal22))
# rtotal24 = np.hstack((r24data, rtotal23))
#
# data0 = np.vstack((rtotal04, rtotal14))
# data = np.vstack((data0, rtotal24))
# print(data.shape)
# hdr['NAXIS1'] = data.shape[1]
# hdr['NAXIS2'] = data.shape[0]
# hduN = fits.PrimaryHDU(data)
# hduN.header = hdr
# # hduN.header['BUNIT'] = '$10^{20} cm^{-2}$'
# hduN.writeto('./infrared.fits', overwrite=True)
# -----------------------------------------

# Post stack bump algorithm stack-bump的后处理方法
# image_file1 = fits.open('/Users/administrator/PMO/Icarus/DataBase/G140_170/OverallData/C18O/G140_170_L2.fits')[0]
# image_file = fits.open('/Users/administrator/PMO/Dorian/314802_moment.fits')[0]
# df = pd.read_csv('/Users/administrator/PMO/Dorian/314802_Pixel_18.csv')
# B = df.loc[:, 'Pixel_L'].values
# L = df.loc[:, 'Pixel_B'].values
# Cmin = df.loc[:, 'Pixel_C0'].values[0]
# Cmax = df.loc[:, 'Pixel_C1'].values[0]
# # hdr = fits.open('/Users/administrator/PMO/Icarus/DataBase/G140_170/OverallData/C18O/G140_170_L2.fits')[0].header
# # Sundries.CalculateMoment('/Users/administrator/PMO/Icarus/DataBase/G140_170/OverallData/C18O/G140_170_L2.fits', 0, (Sundries.c2v(hdr, Cmin),Sundries.c2v(hdr, Cmax)), './L2moment.fits', ty=False)
# image_file = fits.open('/Users/administrator/PMO/Dorian/L2moment.fits')[0]
# data12 = fits.open('/Users/administrator/PMO/Dorian/314802_moment.fits')[0].data
# image_hdr = image_file.header
# data = image_file.data
# p_start = (image_hdr['NAXIS1'] - image_hdr['CRPIX1']) * image_hdr['CDELT1']
# p_end = (0 - image_hdr['CRPIX1']) * image_hdr['CDELT1']
# v_start = ((image_hdr['NAXIS2'] - image_hdr['CRPIX2']) * image_hdr['CDELT2'])
# v_end = ((0 - image_hdr['CRPIX2']) * image_hdr['CDELT2'])
# p_r = (p_start, p_end)
# v_r = (v_start, v_end)

# patus = np.zeros((image_hdr['NAXIS2'], image_hdr['NAXIS1']), dtype=int)
# count = 0
# for i in range(len(L)):
#     patus[B[i], L[i]] = data[B[i], L[i]]
#     count += 1
# patus = np.where(patus == 0, np.nan, patus)

# fig, ax = plt.subplots()
# # clr = ['red', 'orange', 'yellow', 'green', 'cyan', 'blue', 'purple']
# img = ax.imshow(patus, extent=[p_r[1], p_r[0], v_r[1], v_r[0]], cmap='ocean', origin='lower', vmin=0, vmax=4)
# divider = make_axes_locatable(ax)
# cax1 = divider.append_axes("right", size="5%", pad=0)
# cax1 = fig.colorbar(img, cax=cax1)
# cax1.set_label('K km/s')
# # extended_clr = clr * (1200 // len(clr)) + clr[:1200 % len(clr)]
# x = np.arange(p_r[1], p_r[0])
# y = np.arange(v_r[1], v_r[0])
# np.meshgrid(x, y)

# create a zero array

# patus1 = np.zeros_like(data12)
# indices = np.where(~np.isnan(data12))
# for i in range(len(indices[0])):
#     patus1[indices[0][i], indices[1][i]] = 1
# ax.contour(patus1, levels=[0.1], extent=[p_r[1], p_r[0], v_r[1], v_r[0]], colors='darkorange', alpha=0.3, linewidths=0.5)

# for i in tqdm(range(len(id))):
#     patus = np.zeros((image_hdr['NAXIS2'], image_hdr['NAXIS1']), dtype=int)
#     indices = np.where(data == id[i])
#     P = list(zip(indices[1], indices[2]))
#     index = list(set(tuple(i) for i in P))
#     for j in index:
#         patus[j[0]][j[1]] = 1
#     ax.contour(patus, levels=[1], extent=[p_r[1], p_r[0], v_r[1], v_r[0]], colors=extended_clr[i], alpha=0.4, linewidths=0.3)

# indices = np.where(data == 147)
# P = list(zip(indices[1], indices[2]))
# index = list(set(tuple(i) for i in P))
# print(index)
# ax.contour(patus, levels=[1], extent=[p_r[1], p_r[0], v_r[1], v_r[0]], colors='red', alpha=0.7, linewidths=0.5)
# ax.invert_xaxis()
# ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%d°'))
# ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%d°'))
# ax.set_xlabel("Longitude")
# ax.set_ylabel("Latitude")
# fig.savefig('./18within.png', dpi=500)
# fig.savefig('./18within.pdf')
# fig.show()
