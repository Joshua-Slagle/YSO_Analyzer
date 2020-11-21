from astroquery.simbad import Simbad
from astroquery.irsa import Irsa
import astroquery.ibe as IBE
import astropy.coordinates as coord
import astropy.units as u

import matplotlib.pyplot as plt
import matplotlib.markers as markers
import matplotlib.cm as cm
import numpy as np
import math


### WISE Central Wavelengths and Zero Points
wavelengths_um = {
    "w1" : 3.4, 
    "w2" : 4.6,
    "w3" : 12,
    "w4" : 22,
}

zeropoints_jy = {
    "w1" : 309.54,
    "w2" : 171.79,
    "w3" : 31.676,
    "w4" : 8.3635
}


ra_args = []
dec_args = []
radius_args = 0.001 * u.deg

### READING IN FILE
f = open("coordinates.txt", "r")
frame_args = f.readline().rstrip()

index = 0
for x in f:
    y, z = x.split(" ")
    ra_args.append(y.rstrip())
    dec_args.append(z.rstrip())


### WISE
wise_results = []
for index in range(min(len(ra_args), len(dec_args))):
    wise_results.append(Irsa.query_region(coord.SkyCoord(str(ra_args[index] + " " + dec_args[index]), frame=frame_args), radius=radius_args, catalog='allwise_p3as_psd'))
    index+=1
#print(wise_results)

stellar_objects, bands = (len(dec_args), 4)
wmpro = [[1 for j in range(bands)] for i in range(stellar_objects)]
for stellar_object in range(stellar_objects):
    wmpro[stellar_object][0] = wise_results[stellar_object]["w1mpro"]
    wmpro[stellar_object][1] = wise_results[stellar_object]["w2mpro"]
    wmpro[stellar_object][2] = wise_results[stellar_object]["w3mpro"]
    wmpro[stellar_object][3] = wise_results[stellar_object]["w4mpro"]

Fv = [[1 for j in range(bands)] for i in range(stellar_objects)]
for stellar_object in range(stellar_objects):
    Fv[stellar_object][0] = zeropoints_jy["w1"] * pow(10,(-0.4 * wmpro[stellar_object][0]))
    Fv[stellar_object][1] = zeropoints_jy["w2"] * pow(10,(-0.4 * wmpro[stellar_object][1]))
    Fv[stellar_object][2] = zeropoints_jy["w3"] * pow(10,(-0.4 * wmpro[stellar_object][2]))
    Fv[stellar_object][3] = zeropoints_jy["w4"] * pow(10,(-0.4 * wmpro[stellar_object][3]))

### Conversions
jy_to_cgs = pow(10, -23)
speed_of_light_cgs = 2.99792458 * pow(10, 10)
um_to_cm = pow(10, -4)

### Calculating Flux Density for Wavelength
Fy = [[1 for j in range(bands)] for i in range(stellar_objects)]
for stellar_object in range(stellar_objects):
    Fy[stellar_object][0] = Fv[stellar_object][0] * jy_to_cgs * ((speed_of_light_cgs)/pow(um_to_cm * wavelengths_um["w1"], 2))
    Fy[stellar_object][1] = Fv[stellar_object][1] * jy_to_cgs * ((speed_of_light_cgs)/pow(um_to_cm * wavelengths_um["w2"], 2))

    Fy[stellar_object][2] = Fv[stellar_object][2] * jy_to_cgs * ((speed_of_light_cgs)/pow(um_to_cm * wavelengths_um["w3"], 2))

    Fy[stellar_object][3] = Fv[stellar_object][3] * jy_to_cgs * ((speed_of_light_cgs)/pow(um_to_cm * wavelengths_um["w4"], 2))

### Calculating Wavelength times Flux Density for Wavelength
yFy = [[1 for j in range(bands)] for i in range(stellar_objects)]
for stellar_object in range(stellar_objects):
    yFy[stellar_object][0] = Fy[stellar_object][0] * um_to_cm * wavelengths_um["w1"]
    yFy[stellar_object][1] = Fy[stellar_object][1] * um_to_cm * wavelengths_um["w2"]
    yFy[stellar_object][2] = Fy[stellar_object][2] * um_to_cm * wavelengths_um["w3"]
    yFy[stellar_object][3] = Fy[stellar_object][3] * um_to_cm * wavelengths_um["w4"]

### PLotting SEDs
## Chart Details
# Axis Titles
plt.title("log(yFy) vs log(y)")
plt.xlabel("log(y)")
plt.ylabel("log(yFy)")

# Color
colors = cm.rainbow(np.linspace(0, 1, stellar_objects))


## X Data
x = []
x.append(math.log(um_to_cm * wavelengths_um["w1"], 10))
x.append(math.log(um_to_cm * wavelengths_um["w2"], 10))
x.append(math.log(um_to_cm * wavelengths_um["w3"], 10))
x.append(math.log(um_to_cm * wavelengths_um["w4"], 10))

## Y Data
y = [[1 for j in range(bands)] for i in range(stellar_objects)]
for stellar_object in range(stellar_objects):
    y[stellar_object][0] = math.log(yFy[stellar_object][0], 10)
    y[stellar_object][1] = math.log(yFy[stellar_object][1], 10)
    y[stellar_object][2] = math.log(yFy[stellar_object][2], 10)
    y[stellar_object][3] = math.log(yFy[stellar_object][3], 10)
    plt.plot(x, y[stellar_object], color=colors[stellar_object], label=wise_results[stellar_object]["designation"][0].decode("utf-8"), marker="o")

## Plot
#fig, (ax1, ax2) = plt.subplots(1, 2)
plt.legend(loc='best')
plt.tight_layout()
#plt.show()


### CLASSIFICATION
## Denominator
d = x[3] - x[0]

a = []
for stellar_object in range(stellar_objects):
    a.append(y[stellar_object][3] - y[stellar_object][0])

classification = []
for stellar_object in range(stellar_objects):
    alpha = a[stellar_object]/d
    if(alpha > 0.3):
        classification.append("Class I")
    elif (alpha > -0.3):
        classification.append("Flat Spectrum")
    elif (alpha > -1.6):
        classification.append("Class II")
    elif (-1.6 > alpha):
        classification.append("Class III")
    else:
        classification.append("NOT CLASSIFIED")


for stellar_object in range(stellar_objects):
    print(wise_results[stellar_object]["designation"][0].decode("utf-8") + " : " + classification[stellar_object])

### Transition Disk Sources


### Plot at the end
plt.show()


