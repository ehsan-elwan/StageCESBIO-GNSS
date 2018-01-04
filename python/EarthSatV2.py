#!/usr/bin/python

"""
Created on Tue, 09 Jul 2017 11:01:31

Short Description:
::todo
@author: ELWAN

"""
import numpy as np
import ephem
from Download_v2 import *
import openpyxl


def earthsat_builder(df, ob_lat, ob_lon, ob_alt, ddate, cn):
    """

    :param df:
    :param ob_lat:
    :param ob_lon:
    :param ob_alt:
    :param ddate:
    :param cn:
    :return: DataFrame



    Observer: Obs
    - Computes the position of the Body.

    - Uses the date of the observer.

    - Uses the epoch of the observer.

         lon,lat : Latitude and longitude (String)
         elevation : (Int)
         ddate : (String) (date format : '1984/5/30 16:22:56' )
         earthSat : (ephem.Body)

    """
    g = 6.67408e-11
    m_earth = 5.97219e24
    column_to_add = ["el", "az", "rx_lon", "rx_lat", "rx_alt", "cn0", "ddate"]
    for col in range(len(column_to_add)):
        df[column_to_add[col]] = np.nan

    res=pd.DataFrame()
    for rx_lat, rx_lon, rx_alt, rx_date, signal_cn in zip(ob_lat, ob_lon, ob_alt, ddate, cn):
        act_df = df
        ob_lon = np.deg2rad(rx_lon)
        ob_lat = np.deg2rad(rx_lat)
        obs = ephem.Observer()
        obs.lat = np.deg2rad(ob_lat)
        obs.long = np.deg2rad(ob_lon)
        obs.elevation = rx_alt
        for index, row in act_df.iterrows():
            esat = ephem.EarthSatellite()

            esat._e = act_df.at[index, 'Eccentricity']

            esat._M = np.rad2deg(act_df.at[index, 'Mean Anom(rad)'])

            esat._ap = np.rad2deg(act_df.at[index, 'Argument of Perigee(rad)'])

            esat._raan = np.rad2deg(act_df.at[index, 'Right Ascen at Week(rad)'])

            esat._n = np.sqrt(g * m_earth / (act_df.at[index, 'SQRT(A)  (m 1/2)'] ** 2) ** 3) * (86400 / (2 * np.pi))

            esat._epoch = ephem.Date(rx_date)

            esat._inc = np.rad2deg(act_df.at[index, 'Orbital Inclination(rad)'])

            obs.date = ephem.Date(rx_date.strftime('%Y/%-m/%-d %H:%M:%S'))

            esat.compute(obs)

            values = [np.rad2deg(esat.alt), np.rad2deg(esat.az), rx_lon,
                      rx_lat, rx_alt, signal_cn, rx_date]
            for i in range(len(values)):
                act_df.set_value(index, column_to_add[i], values[i])

        res=pd.concat([res, act_df])

    res=res.reset_index(drop=True)

    return df_filter_info(res)


def df_filter_info(df):
    # type: (pd.DataFrame) -> pd.DataFrame
    """
        :param df:
        :return: a filtered DataFrame containing the (necessary data)*
        to calculate specular point coordinates.
        _________________________________________________________________
        necessary data: [RX_lon, RX_lat, RX_alt, dtime, PRN, el, az, cn0]
        -----------------------------------------------------------------

        Take the original data frame and remove the unnecessary data.
    """

    df1 = pd.DataFrame()
    df1['rx_lat'] = df['rx_lat']
    df1['rx_lon'] = df['rx_lon']
    df1['rx_alt'] = df['rx_alt']
    df1['ddate'] = df['ddate']
    df1['prn_id'] = df['prn_id']
    df1['az'] = df['az']
    df1['el'] = df['el']
    df1['cn0'] = df['cn0']

    return df1


if __name__ == "__main__":
    test = OrbitInfo(date=datetime.datetime(2016, 9, 28, 20, 30, 55, 782000))
    # datetime.datetime(2013, 9, 28, 20, 30, 55, 782000)
    df = pd.DataFrame(test.df)
    df2 = earthsat_builder(df, [44.695469,42.695469,42.12345], [-0.854618,-0.954618,0.9123], [1000,3000,1500],
                           [datetime.datetime(2016, 9, 18, 13, 11, 55, 00000),
                            datetime.datetime(2016, 9, 18, 14, 30, 55, 00000),
                            datetime.datetime(2016, 9, 18, 15, 40, 55, 00000)], [np.nan,np.nan,np.nan])
    print(df2)

    #import openpyxl
    #writer = pd.ExcelWriter('output.xlsx')
    #df2.to_excel(writer, 'Sheet1')
    #writer.save()
