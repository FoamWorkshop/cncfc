import numpy as np

def angle(v, v_ref, norm = True):
    ''' normalized angle (0..2pi) between 2 vectors\n
    input: v - ref. vector, v_ref - second vector:\n
         v2
         /
        / angl
       -------> v1
    '''
    angl = np.arctan2(v_ref[1], v_ref[0]) - np.arctan2(v[1], v[0])

    if norm:
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

def normv(v):
    return v/np.linalg.norm(v)

def proj_vector2plane(u, n, y_dir = np.array([0,0,1])):
    '''vector to plane projection.\n
    input: u - vector, n - plane normal, O - attachement point\n
    ref. https://www.maplesoft.com/support/help/maple/view.aspx?path=MathApps%2FProjectionOfVectorOntoPlane'''

    nn = normv(n)
    l_x = normv(np.cross(y_dir, nn))
    l_y = normv(np.cross(l_x, nn))
    u_proj = u - np.dot(u,nn) / np.linalg.norm(nn)**2 * nn

    # print('proj x',)
    # print('proj y',)

    # u_proj /= np.linalg.norm(u_proj)
    return normv(np.hstack((-np.dot(u_proj, l_x), -np.dot(u_proj, l_y))))
    #_p2-u_proj_p1

def transform(spar_data0):
    ''' resultant data collection:\n
    angle, distance, slope vector
    '''
    v_ref = np.array([1,0,0])

    spar_data1 = np.roll(spar_data0, -1, axis=0)
    # print('s input')
    # print(spar_data0)
    # print(spar_data1)
    # print('e input')
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
            PS, slope_v, R = cross_point(P1, P2, Q1, Q2)
            v=proj_vector2plane(Q2-Q1, slope_v, np.array([0, 0, 1]))

            buff_PS.append(PS)
            buff_R.append(R)
            buff_v.append(v)
            buff_slope_v.append(slope_v)

        buff_2PS.append(np.vstack(buff_PS))
        buff_2R.append(np.vstack(buff_R))
        buff_2v.append(np.vstack(buff_v))
        buff_slope_2v.append(np.vstack(buff_slope_v))
#
#angle calc section
#
    np.set_printoptions(threshold='nan')
    a=np.array(buff_2PS)+np.array(buff_slope_2v)
    b=np.roll(a,-1,axis=0)
    # print(a)
    rot_ang=np.zeros(a.shape[:2])
#-->calculate offset for the first idx_row (along spars)
    a_first=np.flipud(a[0])
    b_first=np.roll(a_first,-1, axis=0)

    rot_first_ang=np.zeros(a_first.shape[:1])
    # print(a_first, b_first)
    for i, (var_a, var_b) in enumerate(zip(a_first,b_first)):
        rot_first_ang[i]= angle(var_a, var_b, norm=False)

    rot_first_ang=np.cumsum(np.insert(rot_first_ang[:-1],0,0))
#--<
#-->calculate offsets between section points
    for i, (var_a, var_b) in enumerate(zip(a,b)):
        # print('break')
        for j, (var2a, var2b) in enumerate(zip(var_a, var_b)):
            rot_ang[i,j]= angle(var2a, var2b)

#--<
    # print('rot ang',np.degrees(rot_ang))
    angle_arr = rot_ang
    angle_arr = np.cumsum(np.rot90(rot_ang) , axis=1)
    angle_arr+= np.vstack(rot_first_ang)
    buff_2ang = np.degrees(angle_arr)
#--<

#return ang, - rotation angle
#R,  - surface radius
#PS, - attachement point to axis [0,0,z]
#v - projected slope

    return (buff_2ang,
            np.roll(np.rot90(np.array(buff_2R  )),-1, axis=1),
            np.roll(np.rot90(np.array(buff_2PS )),-1, axis=1),
            np.roll(np.rot90(np.array(buff_2v )),-1, axis=1))


if __name__ == '__main__':

    # cross_point(P1, P2, Q1, Q2)
test_data = np.array([[[-1, -4, 0],[-1, -4, 1], [-2, -2, 2], [-1, -1, 3]],
                      [[4, -1, 0],[4, -1, 1], [2, -1, 2], [1, 0, 3]],
                      [[-2, 2, 0],[-2, 2, 1], [0, 3, 2], [0, 1, 3]]])
# test_data = np.array([[[-2, -2, 2], [-1, -1, 3]],
#                       [[2, -1, 2], [1, 0, 3]],
#                       [[0, 3, 2], [0, 1, 3]]])
# print(test_data)

buff_2ang, buff_2R, buff_2PS, buff_2v = transform(test_data)

print('ang', buff_2ang)
print('buf2', buff_2R[:,:,0])
print('z', buff_2PS[:,:,2])
print('2v', buff_2v  )
# tp, tq = cross_pointvv(P2-P1, Q2-Q1)
# print(P1-tp, Q1-tq) 4.71238898  4.95736764  5.09289536]
# [[ 4.71238898  4.81205763  4.89224248]
#  [ 4.95736764  5.01088791  5.05541292]
#  [ 5.09289536  5.12479942  5.15223156]]

# print(np.degrees(angle(np.array([1,1,0]), np.array([0,0,0]), np.array([-1,-0.9,0]))))
# print(proj_vector2plane(P3, P1, Q2))
# print(proj_vector2plane(np.array([1, -1, 0]), np.array([1, 1, 0]), np.array([0, 0, 1])))
