import pygmt as pygmt
import pandas as pd
#import numpy as np

region_sca = [-20, 40, 55, 73]
##region_sca = [-20, 40, 55, 70]
##region_sca = [-180, 180, -90, 90]

#data_file = pd.read_csv("geoid_calc_scandinavia_GGM03S.csv")
#data_file = pd.read_csv("geoid_calc_scandinavia_whole_Norway_GGM03S.csv")

#data_file = pd.read_csv("GNSS_Norway_geoidheight_GGM03S.csv")
#data_file = pd.read_csv("GNSS_Norway_all_points.csv")

data_file = pd.read_csv("GNSS_Norway_geoidheight_EGM2008_skip1.csv")
latitudes = data_file["latitude"]
longitudes = data_file["longitude"]
elevations = data_file["geoidheight"]
#####

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
    projection="M10i",
)
fig.coast(
    shorelines="1/0.5p",
    region=region_sca,
    frame="ag",
)

fig.colorbar(
    frame=["a5", "x+lGeoidheight", "y+lm"])
fig.show()
##region=[-180, 180, -90, 90]

#Cyl_stere/30/-20/12c