import pandas as panda
import numpy as np

geoid_model_norway = panda.read_csv("./TBA4251_Geoid_Height_Spherical_Harmonic/data/norway_geoid_height.txt", delim_whitespace=True, usecols=["Bredde", "Lengde", "H-orto", "H-ell", "Geoidehøgde"], index_col=False)
geoid_model_norway_Bredde = geoid_model_norway["Bredde"]
geoid_model_norway_Lengde = geoid_model_norway["Lengde"]
geoid_model_norway_H_orto = np.array(geoid_model_norway["H-orto"])
geoid_model_norway_H_ell = np.array(geoid_model_norway["H-ell"])
geoid_model_norway_Geoidehoyde = np.array(geoid_model_norway["Geoidehøgde"])

geoid_model_norway_multiindex = panda.MultiIndex.from_arrays([geoid_model_norway_Bredde, geoid_model_norway_Lengde], names=["Breddegrad", "Lengdegrad"])
geoid_model_norway_dataframe = panda.DataFrame(np.transpose(np.array([geoid_model_norway_H_orto, geoid_model_norway_H_ell, geoid_model_norway_Geoidehoyde])), index=geoid_model_norway_multiindex)


#print(geoid_model_norway_dataframe.loc[61.6929259311394]["H-orto"])
#print(geoid_model_norway_Bredde[0])
#print(geoid_model_norway_H_orto[0])

def test123():
    for i in range(3):
        print(f"For breddegrad: + {geoid_model_norway_Bredde[i]} har vi denne geoidehøyden: {geoid_model_norway_H_ell[i] - geoid_model_norway_H_orto[i]}")

test123()