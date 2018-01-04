#!/usr/bin/python
# coding=utf-8
import cgi
import cgitb

print "Content-type:text/html\r\n\r\n"
print '<meta http-equiv="Content-type" content="text/html;charset=UTF-8">'
from Calc_specular_point import *

cgitb.enable()
data = cgi.FieldStorage()
rx_lat = float(data["rx_lat"].value)
rx_lon = float(data["rx_lon"].value)
rx_alt = float(data["rx_alt"].value)
ddate = data["rx_date"].value
ddate = datetime.datetime.strptime(ddate.replace("/", "-"), "%Y-%m-%d %H:%M")

test = OrbitInfo(ddate)  # Yuma file
df = pd.DataFrame(test.df)  # build a Dataframe
df_spec = earthsat_builder(df, [rx_lat], [rx_lon], [rx_alt],
                           [ddate],[np.nan])  # calc az, el


spec = SpecCalculator(df_spec)
df_spec = spec.spec_calc(10)

