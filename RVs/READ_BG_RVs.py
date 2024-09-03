import numpy as np
from astropy.table import Table
from astropy.time import Time


def get_rolling_means(t,rv,rv_err, dt):
    """time (JD), rv (km/s), dt (days)"""

    bins = np.linspace(t[0],t[-1],int(np.ceil((t[-1]-t[0])/dt)) )
    xpm = (bins[1:]+bins[:-1])/2
    ypm = np.zeros_like(xpm)
    epm = np.zeros_like(xpm)
    npm = np.zeros_like(xpm,dtype=int)

    for i in range(len(bins)-1):
        sel = (t >= bins[i]) & (t < bins[i+1])
        ypm[i] = np.nanmean(rv[sel])
        npm[i] = len(rv[sel])
        
        if npm[i]>0:
            rverr_mean = np.sqrt( np.sum( rv_err[sel]**2 )/len(rv[sel]) )
            rv_std = np.sqrt( np.var(rv[sel]) )
            epm[i] = np.where(rv_std>rverr_mean, rv_std, rverr_mean)
        else:
            epm[i] = 0.0
        

    sel = np.isfinite(ypm)
    return xpm[sel],ypm[sel],epm[sel], npm[sel]

# Lick Observatory
RV_lick = Table.read('RV_Lick_Obs_Vol16.dat',format='ascii.csv',header_start=5,names=['JD','RV'])
RV_lick['RV_err'] = 1.0
RV_lick['source'] = 'Lick v16'
print(RV_lick)

# Cape observatory
RV_cape = Table.read('RV_CapeObs_Vol10.dat',format='ascii.csv',header_start=1,names=['localdate','RV'])
RV_cape['JD']=Time( RV_cape['localdate'], format='iso').jd
RV_cape['RV_err'] = 1.0
RV_cape['source'] = 'Cape v10'
print(RV_cape)

# Bottlinger1911 (Potsdam)
RV_Bot = Table.read('RV_Bottlinger1911.dat',format='ascii.csv',header_start=1,names=['date','RV'])
RV_Bot['JD']=Time(RV_Bot['date'], format='iso').jd
RV_Bot['RV_err'] = 1.0
RV_Bot['source'] = 'Bottlinger1911'
print(RV_Bot)



# Mt Wilson
RV_mw = Table.read('RV_MtWilson_Sanford1933.dat',format='ascii.csv',header_start=2,names=['JD','RV'],)
RV_mw['RV_err'] = 1.0
RV_mw['source'] = 'Sanford1933'
print(RV_mw)

# Adams 1956
RV_a = Table.read('RV_Adams1956.dat',format='ascii.csv',header_start=2,names=['JD','RV'])
RV_a['RV_err'] = 1.0
RV_a['source'] = 'Adams1956'
print(RV_a)

# Weymann
RV_w = Table.read('RV_Weymann1962.dat',format='ascii',names=['JD','RV'])
RV_w['RV_err'] = 0.6
RV_w['source'] = 'Weymann1962'
print(RV_w)

# Boesgaard 1975, 1979
RV_b75 = Table.read('RV_Boesgaard1975.dat',format='ascii.csv',header_start=1,names=['date','RV'])
RV_b75['JD'] = Time( RV_b75['date'],format='iso' ).jd
RV_b75['RV_err'] = 0.6
RV_b75['source'] = 'Boesgaard1975'
print(RV_b75)

RV_b79 = Table.read('RV_Boesgaard1979.dat',format='ascii.csv',header_start=1,names=['JD','RV','RV_err'])
RV_b79['source'] = 'Boesgaard1979'
print(RV_b79)

# Goldberg 1984
RV_g = Table.read('RV_Goldberg_KPNO1984.dat',format='ascii.csv',data_start=2,names=['JD','RV'])
RV_g['RV_err'] = 0.6
RV_g['source'] = 'Goldberg1984'
print(RV_g)

# Smith 1989, 1995
RV_ms = Table.read('RV_1995_AD.dat',format='ascii',names=['HJD','RV','?','residRV'])
RV_ms['JD'] = RV_ms['HJD']+2440000.
RV_ms['RV_err'] = 0.6
RV_ms['source'] = 'Smith1989'
print(RV_ms)

# Oak Ridge 2005
RV_or = Table.read("RV_OakRidge2005.dat",format='ascii',names=['MJD','RV','?'])
RV_or['JD'] = RV_or['MJD']+2400000.
RV_or['RV_err'] = 0.1
RV_or['source'] = 'OakRidge2005'
print(RV_or)

# Gray 2008
RV_gray = Table.read("RV_gray.csv",format='ascii')
RV_gray['RV'] = RV_gray['RVoff']+22.25
RV_gray['RV_err'] = 0.5
RV_gray['source'] = 'Gray2008'
print(RV_gray)

# STELLA 
RV_s = Table.read('RV_STELLA_sorted_2023.dat',format='ascii',names=['JD','RV','RV_err'])
RV_s['source'] = 'STELLA'
print(RV_s)

### COMBINE 
def join_col(col):
    return np.concatenate( (RV_lick[col],
                     RV_cape[col],
                     RV_Bot[col],
                     RV_mw[col],
                     RV_a[col],
                     RV_w[col],
                     RV_b75[col],
                     RV_b79[col],
                     RV_g[col],
                     RV_ms[col],
                     RV_or[col],
                     RV_gray[col],
                     RV_s[col]) )
_t  = join_col('JD')
_rv = join_col('RV')
_rve= join_col('RV_err')
_s  = join_col('source')


# RVall table
RVall = Table([_t,_rv,_rve,_s],names=['JD','RV','RV_err','source'])
RVall.sort('JD')
print("CREATING JOINED TABLE:")
print(RVall)


# Time-binned tables
_t,_rv,_erv,_np=get_rolling_means(RVall['JD'],RVall['RV'],RVall['RV_err'], dt=365.25)
RVarm = Table([_t,_rv,_erv,_np],names=['JD','RV','RV_err','NP'])

_t,_rv,_erv,_np=get_rolling_means(RVall['JD'],RVall['RV'],RVall['RV_err'],dt=10)
RV10 = Table([_t,_rv,_erv,_np],names=['JD','RV','RV_err','NP'])

_t,_rv,_erv,_np=get_rolling_means(RVall['JD'],RVall['RV'],RVall['RV_err'],dt=30)
RV30 = Table([_t,_rv,_erv,_np],names=['JD','RV','RV_err','NP'])

RVall.write('BG_rvall.csv',overwrite=True,format='ascii.csv')
RVarm.write('BG_rv365.csv',overwrite=True,format='ascii.csv')
RV10.write('BG_rv10.csv',overwrite=True,format='ascii.csv')
RV30.write('BG_rv30.csv',overwrite=True,format='ascii.csv')


