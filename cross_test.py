import numpy as np

def angle(v, v_ref):
    ''' normalized angle (0..2pi) between 2 vectors\n
    input: v - ref. vector, v_ref - second vector:\n
         v2
         /
        / angl
       -------> v1
    '''
    angl = np.arctan2(v_ref[1], v_ref[0]) - np.arctan2(v[1], v[0])
    if angl < 0:
        angl += 2 * np.pi

    return angl

def cross_point(P1, P2, Q1, Q2):
    '''attachement point, vector and distance of a segment /z/ perpenticular to 2 vectors\n
    input: P1, P2 - u vector points, Q1, Q2 - v vector points\n
    function returns a tupple: PS attachment point on p line, PSQS vector, length of the segment /z/
    ref. geomalgorithms.com #distance between 3d lines & segments#
    '''
    u=P2-P1
    v=Q2-Q1
    u=u/np.linalg.norm(u)
    v=v/np.linalg.norm(v)
    w0=P1-Q1
    a=np.dot(u,u)
    b=np.dot(u,v)
    c=np.dot(v,v)
    d=np.dot(u,w0)
    e=np.dot(v,w0)
    sc=((b*e)-(c*d))/((a*c)-b**2)
    tc=((a*e)-(b*d))/((a*c)-b**2)

    PS = P1 + sc * u
    QS = Q1 + tc * v
    z = QS - PS
    d = np.linalg.norm(z)

    return (PS, z, d)

def proj_vector2plane( u, n, O = np.array([0, 0, 0])):
    '''vector to plane projection.\n
    input: u - vector, n - plane normal, O - attachement point\n
    ref. https://www.maplesoft.com/support/help/maple/view.aspx?path=MathApps%2FProjectionOfVectorOntoPlane'''

    u_proj = u - np.dot(u,n) / np.linalg.norm(n)**2 * n
    return u_proj + O

def transform(spar_data):
    '''
    '''
    v_ref = np.array([1,0,0])
    angl = np.apply_along_axis(angle, 0, spar_data, v_ref=v_ref)
    print(angl)
    cross_point(P1, P2, Q1, Q2)

test_data = np.arange(27).reshape(3,3,3)
P1=np.array([0,0,0])
P2=np.array([0,-1,0])
# P3=np.array([5,4,2])
Q1=np.array([-2,-1,3])
Q2=np.array([0,0,3])

# P, z, d = cross_point(P1, P2, Q1, Q2)
# print(P, z, d)
print(transform(test_data))
# tp, tq = cross_pointvv(P2-P1, Q2-Q1)
# print(P1-tp, Q1-tq) 4.71238898  4.95736764  5.09289536]
# [[ 4.71238898  4.81205763  4.89224248]
#  [ 4.95736764  5.01088791  5.05541292]
#  [ 5.09289536  5.12479942  5.15223156]]

# print(np.degrees(angle(np.array([1,1,0]), np.array([0,0,0]), np.array([-1,-0.9,0]))))
# print(proj_vector2plane(P3, P1, Q2))
