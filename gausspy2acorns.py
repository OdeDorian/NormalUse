import numpy as np
data = np.loadtxt('/Users/administrator/gausspyplus/example/decomposition_grs/gpy_decomposed/Program1_13_g+_fit_fin_sf-p2_finalized.dat', skiprows=1)
# x represents the pixel index of longitude
x = np.array(data[:, 0])
# y represents the pixel index of latitude
y = np.array(data[:, 1])
# amp represents the intensity
amp = np.array(data[:, 4])
# amp_err represents the intensity error
amp_err = np.array(data[:, 5])
# v for velocity
v = np.array(data[:, 6])
# v_err for velocity error
v_err = np.array(data[:, 7])
# rms for root mean squared noise level
rms = np.array(data[:, 12])
# v_disper stands for velocity dispersion
v_disper = np.array(data[:, 8])
# v_disper_err for velocity dispersion error
v_disper_err = np.array(data[:, 9])
FWHM = 2 * np.sqrt(2 * np.log(2)) * v_disper
FWHM_err = 2 * np.sqrt(2 * np.log(2)) * v_disper_err
print(">>>>>>>>>>>>>>>>>>>>>>>>>>")
stacked_arr = np.column_stack((x, y, amp, amp_err, v, v_err, FWHM, FWHM_err, rms))
np.savetxt('./Output/unassigned_catalogue.dat', stacked_arr, fmt='%10.3f')
print('        Finished')
print('==========================')