import pygmt as pygmt
import pandas as pd
#import numpy as np

region_sca = [-180, 180, -90, 90]


data_file = pd.read_csv("geoid_calc_whole_word_GGM03S.csv")

latitudes = data_file["latitude"]
longitudes = data_file["longitude"]
elevations = data_file["geoidheight"]
####

grid = pygmt.xyz2grd(
    x=longitudes,
    y=latitudes,
    z=elevations,
    spacing=(1,1),
    region=region_sca,
)
fig = pygmt.Figure()
fig.grdimage(
    grid=grid,
    # specify projection within the first plotting method used for this figure
    projection="Q12c",
)
fig.coast(
    shorelines="1/0.5p",
    region=region_sca,
    frame="afg",
)

fig.colorbar(
    frame=["a15", "x+lGeoidheight", "y+lm"])
fig.show()
