import json
from EarthSatV2 import *
from Calc_specular_point import *
from time import sleep


def geojson_to_list(path):
    cords = [[], [], [], [], []]
    with open(path) as f:
        data = json.load(f)

    for feature in data['features']:

        for lon, lat, alt, ddate in feature['geometry']['coordinates']:
            # print lon, lat,alt,ddate
            cords[0].append(float(lat))
            cords[1].append(float(lon))
            cords[2].append(float(alt))
            cords[3].append(datetime.datetime.strptime(ddate, "%Y-%m-%d %H:%M:%S"))
            cords[4].append(np.nan)

    return cords


path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir)) + "/GeoJson/Obslayer.json"
cord = geojson_to_list(path)
test = OrbitInfo(cord[3][0])
df = pd.DataFrame(test.df)
df2 = earthsat_builder(df, cord[0], cord[1], cord[2],
                       cord[3], cord[4])
spec = SpecCalculator(df2)
df_spec = spec.spec_calc(10)
