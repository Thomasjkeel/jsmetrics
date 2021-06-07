import numpy as np
from seawater import dist

def calc_ZWRI(ua, dlon):
    """
        To quantify the response to Arctic sea ice loss, we define a Zonal Wind Response Index (ZWRI) that is calculated
         as the difference in zonally averaged zonal wind response between 30-39◦N and 54-63◦N averaged over 600 to 150 hPa,
          which correspond to the regions with the largest zonal wind responses in the multi-model mean
        From hackathon (Louis and Pheobe)
    """
    latN = np.where( (ua.latitude <= 63) & (ua.latitude >= 54) )[0]
    latS = np.where( (ua.latitude <= 39) & (ua.latitude >= 30) )[0]
    pressi = np.where( (ua.plev/100 <= 600) & (ua.plev/100 >= 150) )[0]
    # Area-weighted, by the lon extent
    lon_km_N = [dist([ua.latitude[i],ua.latitude[i]],[0,dlon])[0][0] for i in latN]
    lon_km_S = [dist([ua.latitude[i],ua.latitude[i]],[0,dlon])[0][0] for i in latS]
#     print(ua.latitude[latN].shape, np.nanmean(ua.ua_diff[pressi,latN],0).shape, len(lon_km_N) )

    ZWRI = np.nansum( lon_km_S*np.nanmean(ua.ua_diff[pressi,latS],0) )/np.nansum(lon_km_S) - np.nansum( lon_km_N*np.nanmean(ua.ua_diff[pressi,latN],0) )/np.nansum(lon_km_N)
#     ZWRI = np.nansum( lon_km_N*np.nanmean(ua.ua_diff[pressi,latN],0) )/np.nansum(lon_km_N) - np.nansum( lon_km_S*np.nanmean(ua.ua_diff[pressi,latS],0) )/np.nansum(lon_km_S)
    return ZWRI


def calc_ZW(ua, dlon):
    """
        Fig.5a of Barnes2015
        Barnes (& Povlani?) 2015
    """
    latN = np.where( (ua.latitude <= 70) & (ua.latitude >= 30) )[0]
    pressi = np.where( (ua.plev/100 <= 600) & (ua.plev/100 >= 150) )[0]
    # Area-weighted, by the lon extent
    lon_km_N = [dist([ua.latitude[i],ua.latitude[i]],[0,dlon])[0][0] for i in latN]
    ZW = np.nansum( lon_km_N*np.nanmean(ua.ua_diff[pressi,latN],0) )/np.nansum(lon_km_N) 
    return ZW



def calc_SPV(ua, dlon):
    """
        The strength of the stratospheric polar vortex (SPV) is computed as the zonal-mean zonal wind averaged over 54-66◦N at 10 hPa
        Hackathon (Louis and Pheobe) or general definition
    """

    lati = np.where( (ua.latitude <= 66) & (ua.latitude >= 54) )[0]
    pressi = np.where( ua.plev/100 == 10 )[0][0]
    # Area-weighted, by the lon extent
    lon_km = [dist([ua.latitude[i],ua.latitude[i]],[0,dlon])[0][0] for i in lati]
    SPV = np.nansum( lon_km*ua.ua_diff[pressi,lati] ) / np.nansum(lon_km)
    return SPV




