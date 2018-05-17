from scipy.spatial import ConvexHull
from shapely.geometry import Polygon, Point
import numpy as np
import sys

#answer file
f = open(sys.argv[1], 'rb')
tmp = f.read().split()
points = np.zeros((int(tmp[0]), 2), dtype=int)
for val in range(0, int((len(tmp)-1)/2)):
    points[val][0] = int(tmp[2*val+1])
    points[val][1] = int(tmp[2*val+2])

f.close()

#students file
f1 = open(sys.argv[2], 'rb')
tmp = f1.read().split()
ans = np.zeros((int(tmp[0]), 2), dtype=int)
for val in range(0, int((len(tmp)-1)/2)):
    ans[val][0] = int(tmp[2*val+1])
    ans[val][1] = int(tmp[2*val+2])

f1.close()

poly = Polygon(points)

correct = True
for i in range(len(ans)):
    if(poly.exterior.distance(Point(ans[i][0], ans[i][1])) > 0.0):
        correct = False

if (correct):
    print(1)