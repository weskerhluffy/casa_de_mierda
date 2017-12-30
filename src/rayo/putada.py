from shapely.geometry import Polygon
from shapely.geometry import Point

# construct a valid polygon (i.e. not self-intersecting) with a hole
ext = [(0, 0), (0, 2), (2, 2), (2, 0), (0, 0)]
int = [(1, 1), (1, 1.5), (1.5, 1.5), (1.5, 1)]
myPoly = Polygon(ext,[int])
print("el poly {}".format(myPoly))

# construct a point that is within the bounds of the exterior ring but inside the hole
myPoint = Point(1.25, 1.25)
print(myPoint.within(myPoly))  # returns 'False'

# construct a point that is inside the polygon (avoiding the hole)
myOtherPoint = Point(0.5, 0.5)
print(myOtherPoint.within(myPoly))  # returns 'True'
