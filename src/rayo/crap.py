from matplotlib import pyplot
from descartes import PolygonPatch
import math
import json

# XXX: https://gis.stackexchange.com/questions/67478/how-to-create-a-circle-vector-layer-with-12-sectors-with-python-pyqgis
# XXX: https://gis.stackexchange.com/questions/93136/how-to-plot-geo-data-using-matplotlib-python0
# XXX: https://stackoverflow.com/questions/17734587/why-is-set-xlim-not-setting-the-x-limits-in-my-figure

BLUE = '#6699cc'
GRAY = '#999999'
data = []
with open('sector.json') as f:
    for line in f:
        data.append(json.loads(line))
#print("datos {}".format(data))


def plotFeature(coordlist, myplot):
    #create a polygon geojson-like feature
    poly = {"type": "Polygon", "coordinates": coordlist}
#    poly = coordlist
    print("malo si si {}".format(poly))
    patch = PolygonPatch(poly, fc=BLUE, ec=BLUE, alpha=0.5, zorder=2)
    #plot it on the graph
    myplot.add_patch(patch)



fig = pyplot.figure(1, dpi=180)
ax = fig.add_subplot(121)

for feature in data[0]['features']:
    caca=feature["geometry"]["coordinates"]
#    caca["coordinates"]=caca["coordinates"][0]
    print("coord {}".format(caca))
    plotFeature(caca, ax)

ax.set_xlim(xmin=-40, xmax=80)
ax.set_ylim(ymin=40, ymax=150)
pyplot.show()
pyplot.savefig("teije.png")
