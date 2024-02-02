import sys
import time
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
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
from matplotlib.patches import Ellipse
from astropy.table import Table
import pandas as pd
import matplotlib.gridspec as gridspec
from tqdm import tqdm

# Prepare for HII region information
# HII = pd.read_csv('/Users/administrator/PMO/Icarus/DataReduction/DataInput/table_irsa_catalog_search_results.csv')
# HII_Lon = HII.loc[:, 'glon'].values
# HII_Lat = HII.loc[:, 'glat'].values
# HII_R = HII.loc[:, 'radius'].values
# HII_HR = HII.loc[:, 'hii_region'].values
# HII_Dist = HII.loc[:, 'd/istance'].values
# HII_Angle = HII_R / (HII_Dist * 1000)
# HIIRow = {'l': HII_Lon, 'b': HII_Lat, 'a': HII_Angle, 'e': HII_HR, 'd': HII_Dist}
# Mask_HII = np.where(~pd.isna(HIIRow['e']))
# HIIRow = {i: HIIRow[i][Mask_HII] for i in HIIRow}
#
# Prepare for Infrared information
# InF = pd.read_csv('/Users/administrator/Documents/Infrared.csv')
# InF_Lon = InF.loc[:, 'glon'].values
# InF_Lat = InF.loc[:, 'glat'].values
# InF_Dist = InF.loc[:, 'dist'].values
# InF_Angle = InF.loc[:, 'angle'].values / 3600
# InFRow = {'l': InF_Lon, 'b': InF_Lat, 'd': InF_Dist, 'a': InF_Angle}
#
# Prepare for SNR information
# SNR = pd.read_csv('/Users/administrator/Daedalus/Background/SNR/Catalogs/SNR.csv')
# SNR_Name = SNR.loc[:, 'Name'].values
# SNR_Lon = SNR.loc[:, 'Longitude'].values
# SNR_Lat = SNR.loc[:, 'Latitude'].values
# SNR_A1 = SNR.loc[:, 'Size1'].values / 60
# SNR_A2 = SNR.loc[:, 'Size2'].values / 60
# SNRRow = {'N': SNR_Name, 'l': SNR_Lon, 'b': SNR_Lat, 'a1': SNR_A1, 'a2': SNR_A2}
#
# # Prepare for Star Formation Regions information
# SFR = pd.read_csv('/Users/administrator/Documents/StarFormation.csv')
# SFR_Name = SFR.loc[:, 'Name'].values
# SFR_Lon = SFR.loc[:, 'Longitude'].values
# SFR_Lat = SFR.loc[:, 'Latitude'].values
# SFR_T = SFR.loc[:, 'Detection'].values
# SFRRow = {'n': SFR_Name, 'l': SFR_Lon, 'b': SFR_Lat, 't': SFR_T}
#
#
# image_file = fits.open('/Users/administrator/PMO/Icarus/DataReduction/DataOutput/Fits/IcarI.fits')[0]
# image_file1 = fits.open('/Users/administrator/PMO/Icarus/DataReduction/DataOutput/Fits/IcarII.fits')[0]
# image_file2 = fits.open('/Users/administrator/PMO/Icarus/DataReduction/DataOutput/Fits/IcarIII.fits')[0]
# image_file3 = fits.open('/Users/administrator/PMO/Icarus/DataReduction/DataOutput/Fits/IcarIIII.fits')[0]
# image_data = image_file.data
# image_data1 = image_file1.data
# image_data2 = image_file2.data
# image_data3 = image_file3.data
# image_hdr = image_file.header
# p_start = (image_hdr['NAXIS1'] - image_hdr['CRPIX1']) * image_hdr['CDELT1']
# p_end = (0 - image_hdr['CRPIX1']) * image_hdr['CDELT1']
# v_start = ((image_hdr['NAXIS2'] - image_hdr['CRPIX2']) * image_hdr['CDELT2'])
# v_end = ((0 - image_hdr['CRPIX2']) * image_hdr['CDELT2'])
# p_r = (p_start, p_end)
# v_r = (v_start, v_end)

# plt.figure()
#
# plt.subplots_adjust(wspace=0, hspace=0)
# gs = gridspec.GridSpec(2, 3, width_ratios=[1, 1, 0.05])
#
# ax00 = plt.subplot(gs[0, 0])
# ax01 = plt.subplot(gs[0, 1])
# ax10 = plt.subplot(gs[1, 0])
# ax11 = plt.subplot(gs[1, 1])
#
# ax2 = plt.subplot(gs[:, 2])
#
# img0 = ax00.imshow(image_data, extent=[p_r[1], p_r[0], v_r[1], v_r[0]], cmap='jet', origin='lower', vmin=0, vmax=30, aspect='auto')
#
# ax00.xaxis.set_ticklabels([])
# ax00.yaxis.set_ticklabels([])
# ax00.tick_params(axis='x', direction='in', bottom=True, top=True)
# ax00.tick_params(axis='y', direction='in', left=True, right=True)
#
# img1 = ax01.imshow(image_data1, extent=[p_r[1], p_r[0], v_r[1], v_r[0]], cmap='jet', origin='lower', vmin=0, vmax=30, aspect='auto')
# ax01.xaxis.set_ticklabels([])
# ax01.yaxis.set_ticklabels([])
# ax01.tick_params(axis='x', direction='in', bottom=True, top=True)
# ax01.tick_params(axis='y', direction='in', left=True, right=True)
#
# img2 = ax10.imshow(image_data2, extent=[p_r[1], p_r[0], v_r[1], v_r[0]], cmap='jet', origin='lower', vmin=0, vmax=30, aspect='auto')
# # ax10.xaxis.set_ticklabels([])
# # ax10.yaxis.set_ticklabels([])
# ax10.tick_params(axis='x', direction='in', bottom=True, top=True)
# ax10.tick_params(axis='y', direction='in', left=True, right=True)
# ax10.yaxis.set_major_formatter(mticker.FormatStrFormatter('%d°'))
# ax10.xaxis.set_major_formatter(mticker.FormatStrFormatter('%d°'))
# ax10.set_xlabel('Longitude')
# ax10.set_ylabel('Latitude')
#
# img3 = ax11.imshow(image_data3, extent=[p_r[1], p_r[0], v_r[1], v_r[0]], cmap='jet', origin='lower', vmin=0, vmax=30, aspect='auto')
# ax11.xaxis.set_ticklabels([])
# ax11.yaxis.set_ticklabels([])
# ax11.tick_params(axis='x', direction='in', bottom=True, top=True)
# ax11.tick_params(axis='y', direction='in', left=True, right=True)
#
# for i in range(len(HIIRow['e'])):
#     circle = patches.Circle((HIIRow['l'][i], HIIRow['b'][i]), HIIRow['a'][i], fill=False, color='red', linewidth=0.6)
#     ax00.add_patch(circle)
#     circle1 = patches.Circle((HIIRow['l'][i], HIIRow['b'][i]), HIIRow['a'][i], fill=False, color='red', linewidth=0.6)
#     ax01.add_patch(circle1)
#     circle2 = patches.Circle((HIIRow['l'][i], HIIRow['b'][i]), HIIRow['a'][i], fill=False, color='red', linewidth=0.6)
#     ax10.add_patch(circle2)
#     circle3 = patches.Circle((HIIRow['l'][i], HIIRow['b'][i]), HIIRow['a'][i], fill=False, color='red', linewidth=0.6)
#     ax11.add_patch(circle3)
# for i in range(len(InFRow['l'])):
#     ellipse = Ellipse(xy=(InFRow['l'][i], InFRow['b'][i]), width=InFRow['a'][i], height=InFRow['a'][i], angle=0, fill=False, color='green', linewidth=0.6)
#     ax00.add_patch(ellipse)
#     ellipse1 = Ellipse(xy=(InFRow['l'][i], InFRow['b'][i]), width=InFRow['a'][i], height=InFRow['a'][i], angle=0, fill=False, color='green', linewidth=0.6)
#     ax01.add_patch(ellipse1)
#     ellipse2 = Ellipse(xy=(InFRow['l'][i], InFRow['b'][i]), width=InFRow['a'][i], height=InFRow['a'][i], angle=0, fill=False, color='green', linewidth=0.6)
#     ax10.add_patch(ellipse2)
#     ellipse3 = Ellipse(xy=(InFRow['l'][i], InFRow['b'][i]), width=InFRow['a'][i], height=InFRow['a'][i], angle=0, fill=False, color='green', linewidth=0.6)
#     ax11.add_patch(ellipse3)
# for i in range(len(SNRRow['l'])):
#     rect = patches.Rectangle((SNRRow['l'][i], SNRRow['b'][i]), SNRRow['a1'][i], SNRRow['a2'][i], angle=0, fill=False, color='yellow', linewidth=0.6)
#     ax00.add_patch(rect)
#     rect1 = patches.Rectangle((SNRRow['l'][i], SNRRow['b'][i]), SNRRow['a1'][i], SNRRow['a2'][i], angle=0, fill=False, color='yellow', linewidth=0.6)
#     ax01.add_patch(rect1)
#     rect2 = patches.Rectangle((SNRRow['l'][i], SNRRow['b'][i]), SNRRow['a1'][i], SNRRow['a2'][i], angle=0, fill=False, color='yellow', linewidth=0.6)
#     ax10.add_patch(rect2)
#     rect3 = patches.Rectangle((SNRRow['l'][i], SNRRow['b'][i]), SNRRow['a1'][i], SNRRow['a2'][i], angle=0, fill=False, color='yellow', linewidth=0.6)
#     ax11.add_patch(rect3)
# for i in range(len(SFRRow['l'])):
#     circle = patches.Circle((SFRRow['l'][i], SFRRow['b'][i]), 0.5, fill=False, color='orange', linewidth=0.6)
#     ax00.add_patch(circle)
#     circle1 = patches.Circle((SFRRow['l'][i], SFRRow['b'][i]), 0.5, fill=False, color='orange', linewidth=0.6)
#     ax01.add_patch(circle1)
#     circle2 = patches.Circle((SFRRow['l'][i], SFRRow['b'][i]), 0.5, fill=False, color='orange', linewidth=0.6)
#     ax10.add_patch(circle2)
#     circle3 = patches.Circle((SFRRow['l'][i], SFRRow['b'][i]), 0.5, fill=False, color='orange', linewidth=0.6)
#     ax11.add_patch(circle3)
#
# plt.colorbar(img0, cax=ax2)
# plt.savefig('/Users/administrator/Documents/Information.png', dpi=500)
# plt.show()


# bx.annotate('R1', xy=(155, 5.5), xytext=(155, 6))
# bx.annotate('R2', xy=(152.5, 5.5), xytext=(152.5, 6))
# bx.annotate('R3', xy=(153.3, 3.5), xytext=(153.7, 3.5))
# rect = patches.Rectangle((153.4, 4.2), 3,  1.7, angle=0, fill=False)
# rect_1 = patches.Rectangle((148, 1.6), 5.2, 5, angle=0, fill=False)
# bx.add_patch(rect)
# bx.add_patch(rect_1)

# OB = pd.read_csv('/Users/administrator/Daedalus/Background/OBstars/LBOBstar.csv')
# OB_Lon = OB.loc[:, 'L'].values
# OB_Lat = OB.loc[:, 'B'].values
# OBRow = {'l': OB_Lon, 'b': OB_Lat}
#
# HII = pd.read_csv('/Users/administrator/Daedalus/Background/SFR & HII/1984HII_Region.csv')
# HII_Name = HII.loc[:, 'ID'].values
# HII_Lon = HII.loc[:, 'glon'].values
# HII_Lat = HII.loc[:, 'glat'].values
# HII_R = HII.loc[:, 'Diam'].values / 60
# HIIRow = {'l': HII_Lon, 'b': HII_Lat, 'r': HII_R, 'N': HII_Name}
#
# image_file = fits.open('/Users/administrator/Daedalus/MC/MWISPData/TwelveData/Moment12_0.fits')[0]
# # image_file1 = fits.open('/Users/administrator/PMO/Icarus/DataReduction/DataOutput/Fits/IcarII.fits')[0]
# # image_file2 = fits.open('/Users/administrator/PMO/Icarus/DataReduction/DataOutput/Fits/IcarIII.fits')[0]
# # image_file3 = fits.open('/Users/administrator/PMO/Icarus/DataReduction/DataOutput/Fits/IcarIIII.fits')[0]
# image_data = image_file.data
# # image_data1 = image_file1.data
# # image_data2 = image_file2.data
# # image_data3 = image_file3.data
# image_hdr = image_file.header
# p_start = (image_hdr['NAXIS1'] - image_hdr['CRPIX1']) * image_hdr['CDELT1']
# p_end = (0 - image_hdr['CRPIX1']) * image_hdr['CDELT1']
# v_start = ((image_hdr['NAXIS2'] - image_hdr['CRPIX2']) * image_hdr['CDELT2'])
# v_end = ((0 - image_hdr['CRPIX2']) * image_hdr['CDELT2'])
# p_r = (p_start, p_end)
# v_r = (v_start, v_end)
#
# fig, ax = plt.subplots()
# img0 = ax.imshow(image_data, extent=[p_r[1], p_r[0], v_r[1], v_r[0]], cmap='gray', origin='lower', vmin=0, vmax=40)
# ax.tick_params(axis='x', direction='in', bottom=True, top=True)
# ax.tick_params(axis='y', direction='in', left=True, right=True)
# for i in range(len(HIIRow['l'])):
#     circle = patches.Circle((HIIRow['l'][i], HIIRow['b'][i]), HIIRow['r'][i], fill=False, color='orangered', linewidth=0.9)
#     ax.annotate(HIIRow['N'][i], xy=(HIIRow['l'][i], HIIRow['b'][i]), size=10, color='saddlebrown')
#     ax.add_patch(circle)
# for i in range(len(SNRRow['l'])):
#     ellipse = patches.Ellipse((SNRRow['l'][i], SNRRow['b'][i]), SNRRow['a1'][i], SNRRow['a1'][i], angle=0, fill=False, color='green', linewidth=0.9)
#     ax.add_patch(ellipse)
#     ax.annotate(SNRRow['N'][i], xy=(SNRRow['l'][i], SNRRow['b'][i]), size=10, color='lime')
# for i in tqdm(range(len(OBRow['l']))):
#     circle = patches.Circle((OBRow['l'][i], OBRow['b'][i]), 0.03, fill=False, color='yellow', linewidth=0.3)
#     # ax.annotate(OBRow['N'][i], xy=(OBRow['l'][i], OBRow['b'][i]), size=10, color='saddlebrown')
#     ax.add_patch(circle)
# ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%d°'))
# ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%d°'))
# ax.set_xlabel('Longitude')
# ax.set_ylabel('Latitude')
# ax.set_title('SNR, HII & OB stars')
# divider = make_axes_locatable(ax)
# cax1 = divider.append_axes("right", size="5%", pad=0)
# cax1 = fig.colorbar(img0, cax=cax1)
# cax1.set_label('K km/s')
# fig.savefig('/Users/administrator/Daedalus/Background/Background.png', dpi=500)
# fig.savefig('/Users/administrator/Daedalus/Background/Background.pdf')
# fig.show()



# 从网页上下载fits数据
data = fits.open('/Users/administrator/Downloads/J_MNRAS_485_2895_catalog.dat.gz.fits')[1].data
# 将数据转换为pandas DataFrame
df = pd.DataFrame(data)
# 将DataFrame保存为CSV文件
df.to_csv('SCOPE_PlanckColdClumps.csv', index=False)

# ra-dec到l-b的转化
# df = pd.read_csv('/Users/administrator/Daedalus/Background/OBstars/OBstarCatalog.csv')
# RA = df.loc[:, 'RAdeg'].values
# Dec = df.loc[:, 'DEdeg'].values
# c = SkyCoord(57.597929, 57.000678, unit='deg')
# l = c.galactic.l.degree
# b = c.galactic.b.degree
# print(l, b)


# 处理WISE图像数据，画图
# 打开FITS文件并获取第一个扩展HDU（主HDU）
# hdu = fits.open('/Users/administrator/PMO/WISE12.fits')[0]
# # 创建一个WCS对象，指定投影方式为右手极坐标系
# wcs = WCS(hdu.header)
# # 创建一个子图，并设置投影为WCS对象
# plt.subplot(projection=wcs)
# plt.imshow(hdu.data, vmin=540, vmax=600, origin="lower", cmap='hot')
# plt.colorbar(pad=0)
# plt.xlabel("Longitude")
# plt.ylabel("Latitude")
# # plt.grid(color="white", ls="solid")
# plt.savefig('./WISE12.png', dpi=500)
# plt.savefig('./WISE12.pdf')
# plt.show()















