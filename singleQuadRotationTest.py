import matplotlib.pyplot as plt
import numpy as np
import math

#My defined functions
def signFuncSohail(a):
  b = 1.0;
  if(a>=0.0): b =  1.0;
  if(a< 0.0): b = -1.0;
  return b;

def InnerProduct(a, b, dim):
  ab = 0.0;
  for i in range(dim):
    ab = ab + a[i]*b[i];
  return ab;

#Allocation 
nodeIDs     = [0, 1, 2, 3];
edgeIDs     = np.array([[0,1], [0,2],  [1,3], [2,3]]);
nodeCoords  = np.array([[0.0,0.0], [1.0,0.0], [0.0,1.0],  [1.0,1.0]]);
edgeTangent = np.array([[0.0,0.0], [0.0,0.0],  [0.0,0.0], [0.0,0.0]]);
edgeCentres = np.array([[0.0,0.0], [0.0,0.0],  [0.0,0.0], [0.0,0.0]]);

nNodes = 4;
nEdges = 4;
ndim   = 2;


#Sohails rotation matrix
#and orthogonal basis vectors 
theta = 20 #Rotation angle degrees

#original basis
e1    = np.array([1.0,0.0]);
e2    = np.array([0.0,1.0]);

#rotated basis
v1    = np.array([0.0,0.0]);
v2    = np.array([0.0,0.0]);
TRotate = np.array([[math.cos(math.radians(theta)), -math.sin(math.radians(theta))]
                  , [math.sin(math.radians(theta)),  math.cos(math.radians(theta))]]);

for i in range(ndim):
  for j in range(ndim):
    v1[i] = v1[i] + TRotate[i][j]*e1[j];
    v2[i] = v2[i] + TRotate[i][j]*e2[j];

#Calculate the edge tangent vectors
for i in range(nEdges):
  for j in range(ndim):
    k = edgeIDs[i][0];
    l = edgeIDs[i][1];
    edgeTangent[i][j] = (nodeCoords[l][j] - nodeCoords[k][j])/5.0; #Scaled so it looks better
    edgeCentres[i][j] = (nodeCoords[l][j] + nodeCoords[k][j])/2.0;

#Flip the edge tangent vectors based on orthogonal basis
for i in range(nEdges):
  S = signFuncSohail(InnerProduct(edgeTangent[i], v1, ndim))*signFuncSohail(InnerProduct(edgeTangent[i], v2, ndim))
  for j in range(ndim):
    edgeTangent[i][j] = S*edgeTangent[i][j]

#Plot the edges
quadShapeID = [0, 1, 3, 2, 0];
x = np.zeros(shape=(nNodes+1));
y = np.zeros(shape=(nNodes+1));
for i in range(nNodes+1):
   x[i] = nodeCoords[quadShapeID[i]][0];
   y[i] = nodeCoords[quadShapeID[i]][1];
plt.plot(x,y);

 
#Plot some arrows indication edge tangents
for i in range(nEdges):
  plt.arrow(edgeCentres[i][0], edgeCentres[i][1], edgeTangent[i][0], edgeTangent[i][1], width = 0.05)

# Showing the graph
plt.show()