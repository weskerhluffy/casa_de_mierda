from shapely.geometry import Polygon
from shapely.geometry import Point
from matplotlib import pyplot
from shapely.geometry import mapping
from descartes.patch import PolygonPatch

BLUE = '#6699cc'
GRAY = '#999999'
# construct a valid polygon (i.e. not self-intersecting) with a hole
ext = [(0, 0), (0, 2), (2, 2), (2, 0), (0, 0)]
inte = [[(1, 1), (1, 1.5), (1.5, 1.5), (1.5, 1), (1, 1)]]
myPoly = Polygon(ext,inte)
print("el poly {}".format(myPoly))

# construct a point that is within the bounds of the exterior ring but inside the hole
myPoint = Point(1.25, 1.25)
print(myPoint.within(myPoly))  # returns 'False'

# construct a point that is inside the polygon (avoiding the hole)
myOtherPoint = Point(0.5, 0.5)
print(myOtherPoint.within(myPoly))  # returns 'True'

fig = pyplot.figure(1, dpi=180)
ax = fig.add_subplot(121)
poly = mapping(myPoly)
patch = PolygonPatch(poly, fc=BLUE, ec=GRAY, alpha=0.5, zorder=2)
ax.add_patch(patch)
ax.set_xlim(xmin=-10, xmax=10)
ax.set_ylim(ymin=-10, ymax=10)
pyplot.show()
