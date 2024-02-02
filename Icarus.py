import numpy as np
import pandas as pd
import os
import sys
import Dorian
from tqdm import tqdm
import warnings


Icarus = Dorian.Dorian()


def Properties(DBScanpath, Gaia):
    print('Making Table  >>>>>>>>>>>')
    warnings.filterwarnings("ignore")
    df = pd.read_csv(DBScanpath)
    LayerSet = df.loc[:, '_idx'].values
    AreaSet = df.loc[:, 'area_exact'].values
    VcenterSet = df.loc[:, 'v_cen'].values
    NameL = df.loc[:, 'x_cen'].values
    NameB = df.loc[:, 'y_cen'].values
    LWth = df.loc[:, 'lineWidth'].values
    PL = df.loc[:, 'peakL'].values
    PB = df.loc[:, 'peakB'].values
    PV = df.loc[:, 'peakV'].values
    df = pd.DataFrame(columns=['Name', 'Area (SquareMin)', 'Vcenter (km/s)', 'Distance (pc)', 'Mass (Solar Mass)', 'NMean (1e20/SquareCM)', 'NMax (1e20/SquareCM)', 'TexMean (K)', 'TexPeak (K)', 'TPLocation(LBV)', 'lineWidth (km/s)'])
    DistanceSet = [np.nan] * len(LayerSet)
    MassSet = [np.nan] * len(LayerSet)
    TexSet = [np.nan] * len(LayerSet)
    CDmeanSet = [np.nan] * len(LayerSet)
    CDmaxSet = [np.nan] * len(LayerSet)
    NameSet = [np.nan] * len(LayerSet)
    TmeanSet = [np.nan] * len(LayerSet)

    for i in tqdm(range(0, len(LayerSet))):
        sys.stdout = open(os.devnull, 'w')
        Icarus.MaskMatch(inpath + '12PartI.fits', inpath + '12PartIdbscanS2P4Con1_Clean.fits', LayerSet[i], outfits + '%i_layer.fits' % LayerSet[i])
        Icarus.CutFits(*np.arange(0, 5), Dimension=3, inpath=outfits + '%i_layer.fits' % LayerSet[i], outpath=outfits + 'SuitFits/%i_suit.fits' % LayerSet[i], Layer=True)
        os.remove(outfits + '%i_layer.fits' % LayerSet[i])
        Tex12 = Icarus.CalculateTex(outfits + 'SuitFits/%i_suit.fits' % LayerSet[i], outfits + '%i_Tex.fits' % LayerSet[i], FITS=False)
        Intensity = Icarus.CalculateMoment(outfits + 'SuitFits/%i_suit.fits' % LayerSet[i], 0, [-90, 30], outfits + 'Moment/%i_Moment.fits' % LayerSet[i], ty=True)
        CD12 = Icarus.Caculate12CD(outfits + 'Moment/%i_Moment.fits' % LayerSet[i], Fits=False)
        Dist = Icarus.Distance(inpath + Gaia, outfits + 'Moment/%i_Moment.fits' % LayerSet[i], outfigures + 'CornerMap/%iCM.pdf' % LayerSet[i], outfigures + 'Distance/%iD.pdf' % LayerSet[i])
        # 需要进行修改
        M12 = Icarus.Mass(Dist, AreaSet[i]/0.25, Intensity)
        DistanceSet[i] = Dist
        MassSet[i] = M12
        TexSet[i] = Tex12[1]
        TmeanSet[i] = Tex12[0]
        CDmeanSet[i] = CD12[0]
        CDmaxSet[i] = CD12[1]
        if NameB[i] > 0:
            NameSet[i] = 'Icarus%.2f+%.2f' % (round(NameL[i], 2), round(NameB[i], 2))
        else:
            NameSet[i] = 'Icarus%.2f%.2f' % (round(NameL[i], 2), round(NameB[i], 2))
        data = pd.Series({
            'Name': NameSet[i],
            'Area (SquareMin)': AreaSet[i],
            'Vcenter (km/s)': round(VcenterSet[i], 3),
            'Distance (pc)': DistanceSet[i],
            'Mass (Solar Mass)': MassSet[i],
            'NMean (1e20/SquareCM)': CDmeanSet[i],
            'NMax (1e20/SquareCM)': CDmaxSet[i],
            'TexMean (K)': TmeanSet[i],
            'TexPeak (K)': TexSet[i],
            'TPLocation(LBV)': '(%i, %i, %i)' % (PL[i], PB[i], PV[i]),
            'lineWidth (km/s)': round(LWth[i], 3)
                })
        df = df._append(data, ignore_index=True)
        df.to_csv(outtable + 'CloudCatalog.csv', index=True)
    sys.stdout = sys.__stdout__
    print('Property Table has been established !\n Job is done !')


current_script_path = os.path.realpath(__file__)
parent_directory = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
outfits = parent_directory + '/Icarus/DataReduction/DataOutput/Fits/'
outfigures = parent_directory + '/Icarus/DataReduction/DataOutput/Figures/'
outtable = parent_directory + '/Icarus/DataReduction/DataOutput/Tables/'
inpath = parent_directory + '/Icarus/DataReduction/DataInput/'
GaiaSet = '1700633435680O-result.csv'

Properties(inpath + 'Book1.csv', GaiaSet)
