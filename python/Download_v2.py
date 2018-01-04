#!/usr/bin/python

"""
Created on Tue, 04 Jul 2017 10:49:26

Short Description:
    OrbitInfo class and functions.

@author: ELWAN

This modules present the object OrbitInfo (constructor, methods)
- Download the orbital plane information of a specific date
- Processing the information's, create a dictionary for every (ID,PRN)
- Create a lite that contains all the previous dictionary's

"""
import urllib2
import pandas as pd
import datetime
import pytz
import os



class OrbitInfo:
    def __init__(self, date):
        self.date = date
        self.path = os.path.dirname(os.path.realpath(__file__))+"/yuma_cache/"
        nnn, yyyy = self.date_to_day_and_year(date)
        self.df = self.downloader(nnn, yyyy)

    def date_to_day_and_year(self, date):
        """ Function to convert Python datetime to a GPS year and day of the year

        Args:
            date (datetime): the needed date

        Returns:
            (int) yyyy: the year
            (int) nnn: the day

        """
        nnn = date.timetuple().tm_yday
        yyyy = date.year

        return nnn, yyyy

    def gps_to_py_time(self, wn, sow, leapsec):
        """ Function to convert GPS WN and sow to a Python dtime

        Args:
            wn (int): GPS week number
            sow (float): GPS Seconds of Week (seconds since last sunday 00:00)
            leapsec (int): Number of leap seconds to be subtracted
                           - https://en.wikipedia.org/wiki/Leap_second
                           - http://www.leapsecond.com/java/gpsclock.htm

        Returns:
            datetime.datetime: Python dtime, in UTC timezone

        """
        epoch = 86400 * (10 * 365 + (1980 - 1969) / 4 + 1 + 6 - 2)
        timestamp = epoch + 86400 * 7 * wn + sow - leapsec
        dtime = pytz.utc.localize(datetime.datetime.utcfromtimestamp(timestamp))

        return dtime

    def yumatocache(self, response, nnn, yyyy):
        f = open(self.path + "yuma_" + str(nnn) + "_" + str(yyyy) + ".txt", "w")
        f.write(str(response))
        f.close()

    def downloader(self, nnn, yyyy):
        """ Function to download and processing the orbital plane information's
        from the "U.S. COAST GUARD NAVIGATION CENTER" website.
        create a dictionary for every (ID,PRN)

        Args:
            nnn (int): Is the number of the almanac from 1 to 999,
            yyyy (int): Is the year starting from 1997 to the current year.
            "If a year is not supplied, it will default to the current year."

        Returns:
            ([dict]) arrweeks: a list of dictionary's

        """
        arrweeks = []
        searched_date = datetime.datetime(yyyy, 1, 1) + datetime.timedelta(nnn - 1)

        try:
            # Try to find a cache file for yuma nnn_yyyy:
            f = open(self.path + "yuma_" + str(nnn) + "_" + str(yyyy) + ".txt", "r")
            print "Found cache file for " + str(searched_date)
        except IOError as e:
            if searched_date >= datetime.datetime.now().replace(hour=0, minute=0, second=0, microsecond=0):
                url = "https://www.navcen.uscg.gov/?pageName=currentAlmanac&format=yuma-txt"
                print("date dans le futur! le fichier Yuma le plus recent ...")
            else:
                print("date : " + str(searched_date))
                url = "http://www.navcen.uscg.gov/?Do=getAlmanac&almanac=" + str(nnn) + \
                      "&year=" + str(yyyy) + "&format=yuma"
            response = urllib2.urlopen(url).read()
            if response == "ERROR: No such file":
                print("Navigation Center website: ERROR: No such file")
            elif "restricted" in response:
                print("blocked by Navigation Center website")
            else:
                print("Caching file: yuma_" + str(nnn) + "_" + str(yyyy) + ".txt")
                # Create cache file for nnn yyyy
                self.yumatocache(response, nnn, yyyy)
        # Read file
        f = open(self.path + "yuma_" + str(nnn) + "_" + str(yyyy) + ".txt", "r")
        flines = f.readlines()
        nb_of_lines = len(flines)
        f.close()
        line_index = 1
        week = int(flines[13].split(':')[1].replace(" ", ""))  # get week number
        sow = float(flines[4].split(':')[1].replace(" ", ""))  # get second of the week
        yuma_valuse = []
        while line_index <= nb_of_lines:
            values = [0, 0, 0, sow, 0, 0, 0, 0, 0, 0, 0, 0, week]
            for index in range(2):
                values[index] = int(flines[line_index + index].split(':')[1].replace(" ", ""))
            for index in range(2, 11):
                values[index] = float(flines[line_index + index].split(':')[1].replace(" ", ""))

            line_index += 15
            yuma_valuse.append(values)
        headers = ['prn_id', 'Health', 'Eccentricity', 'Time of Applicability(s)',
                   'Orbital Inclination(rad)', 'Rate of Right Ascen(r/s)',
                   'SQRT(A)  (m 1/2)', 'Right Ascen at Week(rad)',
                   'Argument of Perigee(rad)', 'Mean Anom(rad)',
                   'Af0(s)', 'Af1(s/s)', 'week']
        df = pd.DataFrame(yuma_valuse, columns=headers)

        return df


if __name__ == "__main__":
    test = OrbitInfo(datetime.datetime(2016, 9, 25, 20, 30, 55))
    print test.df
