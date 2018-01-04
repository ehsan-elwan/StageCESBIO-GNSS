#!/usr/bin/python
# coding=utf-8
"""
Created on Tue, 09 Jul 2017 11:01:31

Short Description:
calculate specular point's, create json files (Ellipse's, point's, and observer)

    We used (pyNMEA2 to analyse the NMEA frame's)

@author: ELWAN

"""

import geojson
import pyproj
import srtm
import os
from geopandas import GeoDataFrame
from shapely.geometry import Point
import numpy as np
from Download_v2 import OrbitInfo
from EarthSatV2 import earthsat_builder
import pandas as pd
import datetime

geod = pyproj.Geod(ellps='WGS84')


class SpecCalculator:
    def __init__(self, df):
        self.df = df
        self.polys = []
        self.path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))+"/"

    def spec_to_cvs(self):
        self.df.to_csv('output.csv', sep='\t', encoding='utf-8')
        pass

    def dframetogeo(self, df):
        geometry = [Point(xy) for xy in zip(df.ell_lat, df.ell_lon)]
        df = df.drop(['ell_lat', 'ell_lon'], axis=1)
        crs = {'init': 'epsg:4326'}
        geo_df = GeoDataFrame(df, crs=crs, geometry=geometry)
        return geo_df

    def spec_to_json(self):
        geo = pd.DataFrame()
        geo['long'] = self.df['lonspec']
        geo['lat'] = self.df['latspec']
        geo['elev'] = self.df['el']
        geo['prn_id'] = self.df['prn_id']
        geo['az'] = self.df['az']
        geo['description'] = self.df['demspec']
        features = []
        geo.apply(lambda x: features.append(
            geojson.Feature(geometry=geojson.Point((x['long'],
                                                    x['lat'],
                                                    x['elev'])),
                            properties=dict(name=x['prn_id'],
                                            az=x['az'],
                                            description=x['description'])))
                  , axis=1)

        with open(self.path + 'GeoJson/SpecPointslayer.json', 'w') as fp:
            geojson.dump(geojson.FeatureCollection(features), fp, sort_keys=True)
        return self.path + 'GeoJson/SpecPointslayer.json'

    def spec_to_Poly_json(self, ell_points_list):
        features = []
        for x in ell_points_list:
            features.append(
                geojson.Feature(geometry=geojson.Polygon([[tuple(l) for l in x[0]]]),
                                properties=""))

        with open(self.path + 'GeoJson/Ellipselayer.json', 'w') as fp:
            geojson.dump(geojson.FeatureCollection(features), fp, sort_keys=True)

        return self.path + 'GeoJson/Ellipselayer.json'

    def obs_to_json(self):

        geo = pd.DataFrame()
        geo['long'] = self.df['rx_lon']
        geo['lat'] = self.df['rx_lat']
        geo['elev'] = self.df['rx_alt']
        geo['description'] = 'Obs'
        features = []
        geo.apply(lambda x: features.append(
            geojson.Feature(geometry=geojson.Point((x['long'],
                                                    x['lat'],
                                                    x['elev'])),
                            properties=dict(description=unicode(x['description'].decode('utf8')))))
                  , axis=1)

        with open(self.path + 'GeoJson/Obslayer.json', 'w') as fp:
            geojson.dump(geojson.FeatureCollection(features), fp, sort_keys=True)
        return self.path + 'GeoJson/Obslayer.json'

    def lat_lon_array_to_dem(self, lat_arr, lon_arr):
        """ Function to retrieve SRTM dem value from latitude and longitude arrays.

        This function calls the SRTM.py_ python library to automatically retrieve
        Digital Elevation Model values from the `SRTM mission`_ from
        a latitude vector and a longitude vector.

        Downloading and caching of SRTM data is handled by SRTM.py_

        Args:
            lat_arr (1D numpy.ndarray): vector of latitudes
            lon_arr (1D numpy.ndarray): vector of longitudes

        Returns:
            1D numpy.ndarray: vector of corresponding DEM values

        .. _SRTM.py: https://pypi.python.org/pypi/SRTM.py
        .. _`SRTM mission`: http://www2.jpl.nasa.gov/srtm/)

        """
        dem_data = srtm.get_data()
        dem_array = np.empty(0, dtype=np.float)
        for elat, elon in zip(lat_arr, lon_arr):
            dem = dem_data.get_elevation(elat, elon)
            if not dem:
                dem_array = np.append(dem_array, [0])
            else:
                dem_array = np.append(dem_array, dem)
        return dem_array

    def pos_to_spec_coords(self, lat, lon, h_msl, elev, azim):

        """ Returns specular point coordinates from RX location and Tx pointing

        .. todo:: Descibe the method

        Args:
            lat (float): Rx latitude [deg]
            lon (float): Rx longitude [deg]
            h_msl (float): Rx height over mean sea level [m] RX_alt
            elev (float): Satellite elevation relative to Rx [deg]
            azim (float): Satellite azimuth relative to Rx [deg]

        Returns:
            (float, float, float):
                - Latitude of specular point
                - Longitude of specular point
                - DEM altitude at specular point location

        """

        # logger.debug('Computing fwd geodetic transform (1st iter)')
        endlon, endlat, _ = geod.fwd(lon, lat, azim,
                                     h_msl / np.tan(np.deg2rad(elev)))

        # logger.debug('Performing DEM lookup (1st iter)')
        dem = self.lat_lon_array_to_dem([endlat], [endlon])

        # logger.debug('Computing fwd geodetic transform (2nd iter)')
        lon_spec, lat_spec, _ = geod.fwd(lon, lat, azim,
                                         (h_msl - dem) / np.tan(np.deg2rad(elev)))

        # logger.debug('Performing DEM lookup (2nd iter)')
        dem_spec = self.lat_lon_array_to_dem([lat_spec], [lon_spec])

        # logger.debug('Computing fwd geodetic transform (3nd iter)')
        lon_spec, lat_spec, _ = geod.fwd(lon, lat, azim,
                                         (h_msl - dem_spec) / np.tan(np.deg2rad(elev)))

        # logger.debug('Performing DEM lookup (3nd iter)')
        dem_spec = self.lat_lon_array_to_dem([lat_spec], [lon_spec])

        return lat_spec, lon_spec, dem_spec

    def satpos_to_sma(self, h, lamb, elev):
        """ Function to compute ground footprint ellipse parameters.

        This function compute the semi major axis :math:`r_a`
        and semi minor axis :math:`r_b` of the ellipse describing
        the first fresnel zone.

        .. math::

           r_b &= \\sqrt{\\frac{\\lambda h}{\\sin(\\epsilon')} +
           \\left(\\frac{\\lambda}{2 \\sin(\\epsilon')}\\right)^2}

           r_a &= \\frac{b}{\\sin(\\epsilon')}

        with :math:`\\lambda` the wavelength [m],
        :math:`h` the receiver height [m] and
        :math:`\epsilon'` the satellite elevation seen from the specular
        reflection point [deg] (i.e., corresponding to the reflection angle)

        References:


        Args:
            h (float): Receptor height [m]
            lamb (float): Signal wavelength :math:`\\lambda` [m]
            elev (float): Satellite elevation [deg]

        Returns:
            (float, float):
              - ra: ellipse semi-minor axis [m]
              - rb: ellipse semi-major axis [m]

        Notes:
            The approach used in this routine is based on the equations from
            [ROU2014]_.
            According to A. Camps (See poster from IGARSS 2016),
            the actual zone contributing for the reflection in case
            of coherent scattering might be smaller than the usually considered
            first Fresnel zone. In that case, the ellipse parameters should be
            multiplied by 0.56 to have a better description of the active area.

        .. [ROU2014] N. Roussel, F. Frappart, G. Ramillien, J. Darrozes,
           C. Desjardins, P. Gegout, F. Pérosanz, and R. Biancale,
           “Simulations of direct and reflected wave trajectories for ground-based
           GNSS-R experiments,”
           Geosci. Model Dev., vol. 7, no. 5, pp. 2261–2279, Oct. 2014.

        """
        h_clip = np.clip(h, 0, np.inf)  # clip height to avoid negative values
        rb = np.sqrt(lamb * h_clip / np.sin(np.deg2rad(elev)) +
                     np.power(lamb / 2 / np.sin(np.deg2rad(elev)), 2))
        ra = rb / np.sin(np.deg2rad(elev))
        return ra, rb

    def ellipse(self, lat_s, lon_s, h_og, elev, azim, lambd, vertices=10):
        """ Function to compute the ellipse footprint coordinates

        Args:
            lat_s (float): Latitude of the specular point center [deg]
            lon_s (float): Longitude of the specular point center [deg]
            h_og (float): vertical distance between receiver and surface [m]
            elev (float): satellite elevation angle [deg]
            azim (float): azimuth angle [deg]
            lambd (float): Signal wavelength :math:`\\lambda` [m]
            vertices (int): number of segments in the resulting polygon

        Returns:
            numpy.ndarray: array of lat_s-lon coordinate describing the footprint
        """

        s_maj_axis, s_min_axis = self.satpos_to_sma(h_og, lambd, elev)
        rot = np.deg2rad(-azim)
        angle1, angle2, lat_onv = geod.inv(lon_s, lat_s, lon_s, lat_s + 1)
        angle1, angle2, lon_onv = geod.inv(lon_s, lat_s, lon_s + 1, lat_s)
        coordlist = []
        for angle in np.arange(0, 360.01, 360. / vertices):
            y = s_maj_axis * np.cos(np.deg2rad(angle))
            x = s_min_axis * np.sin(np.deg2rad(angle))
            el_lng = (x * np.cos(rot) - y * np.sin(rot)) / lon_onv
            el_lat = (y * np.cos(rot) + x * np.sin(rot)) / lat_onv
            coordlist.append([lon_s + el_lng, lat_s + el_lat])
        return np.asarray(coordlist)

    def spec_calc(self, el_angle):
        """
        Calculate specular point's coordinates and dem from RX location and Tx pointing

        :rtype: DataFrame
        """

        self.df = pd.concat([self.df, pd.DataFrame(columns=['latspec', 'lonspec',
                                                            'demspec'])])
        self.df = self.df[self.df.el > el_angle]

        polys = []
        lambda_ellipse = 3e8 / 1.57542e9
        self.df = self.df.reset_index(drop=True)
        n = self.df.shape[0]
        # for index, row in self.df.iterrows():
        for index in range(n):
            rx_latitude = self.df.at[index, 'rx_lat']
            rx_longitude = self.df.at[index, 'rx_lon']
            rx_altitude = self.df.at[index, 'rx_alt']
            prn_id = self.df.at[index, 'prn_id']
            # ddate = self.df.at[index, 'ddate']
            el = self.df.at[index, 'el']
            az = self.df.at[index, 'az']
            latspec, lonspec, demspec = self.pos_to_spec_coords(rx_latitude, rx_longitude, rx_altitude, el, az)
            self.df.set_value(index, 'latspec', latspec)
            self.df.set_value(index, 'lonspec', lonspec)
            self.df.set_value(index, 'demspec', demspec[0])

            ellipse_point = self.ellipse(latspec, lonspec, rx_altitude - demspec[0], el, az, lambda_ellipse,
                                         vertices=10)

            polys.append([ellipse_point, prn_id])
        self.spec_to_Poly_json(polys)
        self.spec_to_json()
        self.obs_to_json()

        return self.df.set_index(['rx_lat', 'rx_lon', 'rx_alt', 'ddate', 'prn_id'])


if __name__ == "__main__":
    ####### OFFLINE:
    test = OrbitInfo(datetime.datetime(2016, 9, 18, 00, 00, 55, 00000))  # Yuma file
    df = pd.DataFrame(test.df)  # build a Dataframe
    df_spec = earthsat_builder(df, [43.11702412135048, 43.13306116240612, 43.113014204188914, 43.104993581605505,
                                    43.12905229628564,
                                    43.18515250937298, 43.249203966977845, 43.265206318396025, 43.27720532212024,
                                    43.28520334369384,
                                    43.30919109985686, 43.29320031385282, 43.25320494908846, 43.205175817237304,
                                    43.1450861841603,
                                    43.08493742707592, 43.07691312608711, 43.04881979669318, 43.0287452513488,
                                    42.956422511073335,
                                    42.89206418807337, 42.92827401776912],
                               [1.7633056640625, 1.6973876953125, 1.636962890625, 1.5435791015625, 1.4996337890624998,
                                1.4666748046875, 1.4501953125, 1.395263671875, 1.3238525390624998, 1.2249755859375,
                                1.16455078125, 1.065673828125, 0.999755859375, 0.9667968749999999, 0.9667968749999999,
                                1.0491943359375, 1.131591796875, 1.1920166015625, 1.241455078125, 1.2579345703125,
                                1.3897705078125, 1.4996337890624998, ],
                               [1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
                                1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000],
                               [datetime.datetime(2016, 9, 18, 9, 11, 55, 00000),
                                datetime.datetime(2016, 9, 18, 9, 30, 55, 00000),
                                datetime.datetime(2016, 9, 18, 10, 11, 55, 00000),
                                datetime.datetime(2016, 9, 18, 10, 11, 55, 00000),
                                datetime.datetime(2016, 9, 18, 11, 11, 55, 00000),
                                datetime.datetime(2016, 9, 18, 11, 11, 55, 00000),
                                datetime.datetime(2016, 9, 18, 12, 11, 55, 00000),
                                datetime.datetime(2016, 9, 18, 12, 11, 55, 00000),
                                datetime.datetime(2016, 9, 18, 13, 11, 55, 00000),
                                datetime.datetime(2016, 9, 18, 13, 11, 55, 00000),
                                datetime.datetime(2016, 9, 18, 14, 11, 55, 00000),
                                datetime.datetime(2016, 9, 18, 14, 11, 55, 00000),
                                datetime.datetime(2016, 9, 18, 15, 11, 55, 00000),
                                datetime.datetime(2016, 9, 18, 15, 11, 55, 00000),
                                datetime.datetime(2016, 9, 18, 16, 11, 55, 00000),
                                datetime.datetime(2016, 9, 18, 16, 11, 55, 00000),
                                datetime.datetime(2016, 9, 18, 17, 11, 55, 00000),
                                datetime.datetime(2016, 9, 18, 17, 11, 55, 00000),
                                datetime.datetime(2016, 9, 18, 18, 11, 55, 00000),
                                datetime.datetime(2016, 9, 18, 18, 11, 55, 00000),
                                datetime.datetime(2016, 9, 18, 19, 11, 55, 00000),
                                datetime.datetime(2016, 9, 18, 19, 11, 55, 00000)
                                ]

                               , [np.nan, np.nan, np.nan, np.nan, np.nan,
                                  np.nan, np.nan, np.nan, np.nan, np.nan,
                                  np.nan, np.nan, np.nan, np.nan, np.nan,
                                  np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])

    spec = SpecCalculator(df_spec)
    df_spec = spec.spec_calc(10)
    print spec.df
