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
    slope_v = QS - PS
    d = np.linalg.norm(slope_v)
    return (PS, slope_v, d)

def proj_vector2plane( u, n, O = np.array([0, 0, 0])):
    '''vector to plane projection.\n
    input: u - vector, n - plane normal, O - attachement point\n
    ref. https://www.maplesoft.com/support/help/maple/view.aspx?path=MathApps%2FProjectionOfVectorOntoPlane'''

    u_proj = u - np.dot(u,n) / np.linalg.norm(n)**2 * n
    return u_proj + O


def transform(spar_data0):
    ''' resultant data collection:\n
    angle, distance, slope vector
    '''
    v_ref = np.array([1,0,0])

    spar_data1 = np.roll(spar_data0, -1, axis=0)
    print('s input')
    print(spar_data0)
    print(spar_data1)
    print('e input')
    segment_vect = spar_data1 - spar_data0

    rot_P1 = np.zeros_like(spar_data0)
    rot_P2 = np.ones_like(spar_data0) * [0,0,1]
    buff_1 = []
    buff_2PS = []
    buff_2v = []
    buff_2R = []
    buff_2ang = []
    buff_slope_2v =[]
    for (sP1, sP2, sQ1, sQ2) in zip(rot_P1, rot_P2, spar_data0, spar_data1):
        buff_PS = []
        buff_v = []
        buff_R = []
        buff_ang = []
        buff_slope_v = []

        for (P1, P2, Q1, Q2) in zip(sP1, sP2, sQ1, sQ2):
            # print(P1, P2, Q1, Q2)
            PS, slope_v, R = cross_point(P1, P2, Q1, Q2)
            # print('points: ',P1, P2, Q1, Q2)
            # print('result: ',PS, slope_v, R)
            v = proj_vector2plane(Q2-Q1, slope_v,PS)
            buff_PS.append(PS)
            buff_R.append(R)
            buff_v.append(v)
            buff_slope_v.append(slope_v)


        buff_2PS.append(np.vstack(buff_PS))
        buff_2R.append(np.vstack(buff_R))
        buff_2v.append(np.vstack(buff_v))
        buff_slope_2v.append(np.vstack(buff_slope_v))
        # buff_2ang.append(np.vstack(buff_ang))

    np.set_printoptions(threshold='nan')
    a=np.array(buff_2PS)+np.array(buff_slope_2v)
    b=np.roll(a,-1,axis=0)
    print(a)
    rot_ang=np.zeros(a.shape[:2])
    # print('zero',rot_ang)
#-->calculate offset for the first idx_row (along spars)
    a_first=a[0]
    b_first=np.roll(a[0],-1, axis=0)
    rot_first_ang=np.zeros(a_first.shape[:1])
    # print(a_first, b_first)
    for i, (var_a, var_b) in enumerate(zip(a_first,b_first)):
        rot_first_ang[i]= angle(var_b, var_a)
    rot_first_ang=np.insert(np.cumsum(rot_first_ang)[1:],0,0)
#--<

#-->calculate offsets between section points
    for i, (var_a, var_b) in enumerate(zip(a,b)):
        # print('break')
        for j, (var2a, var2b) in enumerate(zip(var_a, var_b)):
            # print(var2a, var2b)
            rot_ang[i,j]= angle(var2a, var2b)
            # print(angle(var2b, var2a))
            # print(var2b)
    angle_arr = np.insert(np.cumsum(np.degrees(rot_ang), axis=0),0,rot_first_ang,axis=0)
    # b=np.rot90(a)
    # c=np.roll(b,-1,axis=0)
    # print(b)
    # print(c)
    # print('->2',np.roll(np.array(buff_2PS+slope_v),1,axis=0))

    # for (PS, slope_v) in zip(buff_2PS, buff_slope_2v):
                    # print(P1, P2, Q1, Q2)
        # arm = PS + slope_v
        # arm_2 = np.roll(arm, -1, axis=0)
        # print(arm, arm_2)
                    # ang = angle(arm, arm_2)
                    # print(np.degrees(ang))
                    # buff_ang.append(ang)


#return ang, - rotation angle
#R,  - surface radius
#PS, - attachement point to axis [0,0,z]
#v - projected slope
    return (buff_2ang, buff_2R, buff_2PS, buff_2v)

    # cross_point(P1, P2, Q1, Q2)
test_data = np.array([[[-1,-1,0], [-1,-1,1]],
                       [[1,0,0], [1,0,1]],
                       [[0,1,0], [0,1,1]]])

# print(test_data)

buff_2ang, buff_2R, buff_2PS, buff_2v = transform(test_data)

# print('ang',np.degrees(buff_2ang))
# buff_2R
# buff_2PS
# buff_2v
# tp, tq = cross_pointvv(P2-P1, Q2-Q1)
# print(P1-tp, Q1-tq) 4.71238898  4.95736764  5.09289536]
# [[ 4.71238898  4.81205763  4.89224248]
#  [ 4.95736764  5.01088791  5.05541292]
#  [ 5.09289536  5.12479942  5.15223156]]

# print(np.degrees(angle(np.array([1,1,0]), np.array([0,0,0]), np.array([-1,-0.9,0]))))
# print(proj_vector2plane(P3, P1, Q2))
