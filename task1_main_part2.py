import pandas as panda
from task1_main_part1 import N_gravemetric_total_sum

geoid_model_norway = panda.read_csv("./TBA4251_Geoid_Height_Spherical_Harmonic/data/norway_geoid_height.txt", delim_whitespace=True, usecols=["Bredde", "Lengde", "H-orto", "H-ell", "Geoidehøgde"], index_col=False)
geoid_model_norway_Bredde = geoid_model_norway["Bredde"]
geoid_model_norway_Lengde = geoid_model_norway["Lengde"]
geoid_model_norway_geoid_height = geoid_model_norway["Geoidehøgde"]

GGM03S_model_NMAX = 180
EGM2008_model_NMAX = 2190

def calc_geoid_for_Norway_GGM03():
    data = ['latitude,longitude,geoidheight']
    latitudes = geoid_model_norway_Bredde
    longitudes = geoid_model_norway_Lengde
    for i in range (0, len(geoid_model_norway_Bredde), 2):
        print(latitudes[i], longitudes[i])
        geoid_height_calculations = N_gravemetric_total_sum(latitudes[i], longitudes[i], GGM03S_model_NMAX)
        data.append(str(latitudes[i]) + ',' + str(longitudes[i]) + ',' + str(geoid_height_calculations))
    with open('GNSS_Norway_geoidheight_GGM03S_skip1.csv', 'w') as new_file:
        new_file.write('\n'.join(data))

#Remember to change the formulas in main1_part1 before using it for EGM2008
def calc_geoid_for_Norway_EGM2008():
    data = ['latitude,longitude,geoidheight']
    latitudes = geoid_model_norway_Bredde
    longitudes = geoid_model_norway_Lengde

    for i in range (0, len(geoid_model_norway_Bredde), 2):
        print(latitudes[i], longitudes[i])
        geoid_height_calculations = N_gravemetric_total_sum(latitudes[i], longitudes[i], EGM2008_model_NMAX)
        data.append(str(latitudes[i]) + ',' + str(longitudes[i]) + ',' + str(geoid_height_calculations))
    with open('GNSS_Norway_geoidheight_EGM2008.csv', 'w') as new_file:
        new_file.write('\n'.join(data))

def get_every_other_value():
    data = ['latitude,longitude,geoidheight']
    latitudes = geoid_model_norway_Bredde
    longitudes = geoid_model_norway_Lengde
    geoidhøgde = geoid_model_norway_geoid_height
    for i in range(len(geoid_model_norway_Bredde)):
        print(latitudes[i], longitudes[i])
        data.append(str(latitudes[i]) + ',' + str(longitudes[i]) + ',' + str(geoidhøgde[i]))
    with open('GNSS_Norway_all_points.csv', 'w') as new_file:
        new_file.write('\n'.join(data))

#calc_geoid_for_Norway_GGM03()
#calc_geoid_for_Norway_EGM2008()
get_every_other_value()