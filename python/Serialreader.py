#!/usr/bin/python
# coding=utf-8

"""
Created on Tue, 16 Jul 2017 14:01:31

Short Description:
Get satellite information's through serial port with ublox receptor :
[azimuth, elevation, position],from NMEA frame's type (GSV, GGA, RMC)
calculate specular point's, create json file's (Ellipse's, point's, and observer)


@author: ELWAN

"""
from os.path import isfile
import cgi, cgitb
import geojson
import numpy as np
import pynmea2
import pyproj
import serial
import os
import datetime
import pandas as pd
import os.path

import srtm


"""    """
form = cgi.FieldStorage()
geod = pyproj.Geod(ellps='WGS84')

lambda_ellipse = 3e8 / 1.57542e9


def lat_lon_array_to_dem(lat_arr, lon_arr):
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


def pos_to_spec_coords(lat, lon, h_msl, elev, azim):
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
    dem = lat_lon_array_to_dem([endlat], [endlon])

    # logger.debug('Computing fwd geodetic transform (2nd iter)')
    lon_spec, lat_spec, _ = geod.fwd(lon, lat, azim,
                                     (h_msl - dem) / np.tan(np.deg2rad(elev)))

    # logger.debug('Performing DEM lookup (2nd iter)')
    dem_spec = lat_lon_array_to_dem([lat_spec], [lon_spec])

    # logger.debug('Computing fwd geodetic transform (3nd iter)')
    lon_spec, lat_spec, _ = geod.fwd(lon, lat, azim,
                                     (h_msl - dem_spec) / np.tan(np.deg2rad(elev)))

    # logger.debug('Performing DEM lookup (3nd iter)')
    dem_spec = lat_lon_array_to_dem([lat_spec], [lon_spec])

    return lat_spec, lon_spec, dem_spec


def satpos_to_sma(h, lamb, elev):
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


def spec_to_Poly_json(ell_points_list):
    features = []
    path = os.path.dirname(os.path.realpath(__file__)) + "/"
    for x in ell_points_list:
        features.append(
            geojson.Feature(geometry=geojson.Polygon([[tuple(l) for l in x[0]]]),
                            properties=""))

    with open(path + 'GeoJson/Ellipselayer.json', 'w') as fp:
        geojson.dump(geojson.FeatureCollection(features), fp, sort_keys=True)

    return path + 'GeoJson/Ellipselayer.json'


def ellipse(lat_s, lon_s, h_og, elev, azim, lambd, vertices=10):
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

    s_maj_axis, s_min_axis = satpos_to_sma(h_og, lambd, elev)
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


def spec_to_json(df):
    geo = pd.DataFrame()
    path = os.path.dirname(os.path.realpath(__file__)) + "/"
    geo['long'] = df['lonspec']
    geo['lat'] = df['latspec']
    geo['elev'] = df['el']
    geo['PRN'] = df['prn_id']
    geo['az'] = df['az']
    geo['description'] = df['demspec']
    features = []
    geo.apply(lambda x: features.append(
        geojson.Feature(geometry=geojson.Point((x['long'],
                                                x['lat'],
                                                x['elev'])),
                        properties=dict(name=x['PRN'],
                                        az=x['az'],
                                        description=x['description'])))
              , axis=1)

    with open(path + 'GeoJson/SpecPointslayer.json', 'w') as fp:
        geojson.dump(geojson.FeatureCollection(features), fp, sort_keys=True)
    return path + 'GeoJson/SpecPointslayer.json'


def obs_to_json(df):
    geo = pd.DataFrame()
    path = os.path.dirname(os.path.realpath(__file__)) + "/"
    geo['long'] = df['rx_lon']
    geo['lat'] = df['rx_lat']
    geo['elev'] = df['rx_alt']
    geo['description'] = 'Obs'
    features = []
    geo.apply(lambda x: features.append(
        geojson.Feature(geometry=geojson.Point((x['long'],
                                                x['lat'],
                                                x['elev'])),
                        properties=dict(description=unicode(x['description'].decode('utf8')))))
              , axis=1)

    with open(path + 'GeoJson/Obslayer.json', 'w') as fp:
        geojson.dump(geojson.FeatureCollection(features), fp, sort_keys=True)
    return path + 'GeoJson/Obslayer.json'


def start_serial(port, baudrate, timeout):
    flag_second = False
    try:
        ser = serial.Serial(
            port=port,
            baudrate=baudrate,
            timeout=timeout
        )
    except serial.serialutil.SerialException as e:
        print e
        ser = None
        exit(2)

    if ser.isOpen():
        try:
            ser.readline()
            path = os.path.dirname(os.path.realpath(__file__)) + "/"  # create cache file
            filename = str(datetime.datetime.now().strftime('%Y/%-m/%-d %H:%M:%S')).replace('/', '-')
            f = open(path + "Ublox_cache/ublox" + filename + ".txt", "w")

            df = pd.DataFrame(columns=["rx_lat", "rx_lon", "rx_alt", "rx_ddate", "prn_id",
                                       "az", "el", "cn0", "latspec", "lonspec", "demspec"])

            info = [np.nan for n in range(df.shape[1])]
            polys = []
            while ser.isOpen():
                ublox_data = ser.readline()
                result = pynmea2.re.search('\$GP(.*)$', ublox_data)
                if hasattr(result, 'group'):
                    frame = str(pynmea2.NMEASentence.parse(result.group(0)))
                    f.write(str(frame) + " \n")

                    if "$GPRMC" in frame:
                        rmc_frame = pynmea2.parse(frame)
                        try:
                            info[0] = float(rmc_frame.latitude)
                        except:
                            pass
                        try:
                            info[1] = float(rmc_frame.longitude)
                        except:
                            pass
                        try:
                            info[3] = datetime.datetime.combine(rmc_frame.datestamp, rmc_frame.timestamp)
                            # print info[3]
                        except:
                            print "erreur"

                    elif "$GPGGA" in frame:
                        gga_frame = pynmea2.parse(frame)
                        try:
                            info[2] = float(gga_frame.altitude)
                        except:
                            pass
                    elif "$GPGSV" in frame:
                        gsv_frame = pynmea2.parse(frame)
                        org_info = info
                        if int(gsv_frame.msg_num) * 4 > gsv_frame.num_sv_in_view:
                            nb_sat = int(gsv_frame.num_sv_in_view) % 4
                        else:
                            nb_sat = 4

                        for j in range(1, nb_sat + 1):  # Nb of msg's
                            try:
                                info[4] = int(getattr(gsv_frame, "sv_prn_num_" + str(j)))
                            except:
                                pass
                            try:
                                info[5] = float((getattr(gsv_frame, "azimuth_" + str(j))))
                            except:
                                pass
                            try:
                                info[6] = float((getattr(gsv_frame, "elevation_deg_" + str(j))))
                            except:
                                pass
                            try:
                                info[7] = float(getattr(gsv_frame, "snr_" + str(j)))
                            except:
                                pass
                            if info[3].microsecond == 0:
                                info[8], info[9], info[10] = pos_to_spec_coords(info[0], info[1], info[2], info[6],
                                                                                info[5])
                                info[10] = info[10][0]
                                # filter elevation angel
                                if info[6] > 20:
                                    df.loc[df.shape[0]] = info

                                    ellipse_point = ellipse(info[8], info[9], info[2] - info[10], info[6], info[5],
                                                            lambda_ellipse,
                                                            vertices=10)
                                    # print ellipse_point
                                    info = org_info
                                    flag_second = True
                                    polys.append([ellipse_point, info[4]])

                    elif "$GPGLL" in frame:
                        # Update/write my geojson file here
                        if flag_second:
                            df.to_csv(path + "Specular_cache/" + "Spec_" + filename + ".csv", mode='a',
                                      header=False)
                            obs_to_json(df)
                            spec_to_json(df)
                            spec_to_Poly_json(polys)

                            polys = []

                            df = pd.DataFrame(columns=["rx_lat", "rx_lon", "rx_alt", "rx_ddate", "prn_id",
                                                       "az", "el", "cn0", "latspec", "lonspec", "demspec"])
                            flag_second = False

                if os.path.isfile(path + "offline"):
                    os.remove(path + "offline")
                    f.close()
                    exit(2)

            f.close()

        except serial.serialutil.SerialException as e:
            print e
            exit(3)


if __name__ == "__main__":
    start_serial('/dev/ttyACM0', 9600, 1)
