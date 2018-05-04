import numpy as np
import ezdxf
import collections
import re
from scipy import spatial
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



class chain():

    def __init__(self, fname_dxf):
        '''arr_sgm - segments array [[[x,y,z],[x1,y1,z1]]...]
           arr_prp - properties vector [l0, l0, ... l1]
           dct_prp - properties dictionary {l0:{...}, l1:{...}'''

        seg_dt = np.dtype([('seg', float, (2,3)), ('prp', int, 1)])

        self.seg_dt = seg_dt
        self.dwg = ezdxf.readfile(fname_dxf)
        self.seg_arr = np.array([], dtype = seg_dt)
        self.prp_dict = collections.OrderedDict({'glob':{}, 'loc':{}})
        self.model = collections.OrderedDict({})

    def SegArr2Poly(self, seg):
        return np.vstack([seg[:,0], seg[-1,1]])

    def flip(self, l , o):
        tmp = l
        if o:
            tmp['seg'] = np.roll( l['seg'], 1, axis = 0)
            return tmp
        else:
            return l

    def OrderSegArr(self, seg):
        idx_arr = np.argsort(np.linalg.norm(seg['seg'], axis = 2), axis=1)
        return np.stack([ self.flip(l, idx[0]) for l, idx in zip(seg, idx_arr)], axis = 0)

    def NearestPt(self, data, ref):
        tree = spatial.KDTree(data)
        return tree.query(ref)

    def NearestSeg(self, seg_arr, ref):
        dA, iA= self.NearestPt(seg_arr['seg'][:,0,:], ref)
        dB, iB= self.NearestPt(seg_arr['seg'][:,1,:], ref)
        if dA<=dB:
            return seg_arr[iA], iA
        else:
            return self.flip(seg_arr[iB], 1), iB

    def SortSegArr(self, buf_data, pt):
        # print('entered sort seg arr')
        res_data = buf_data
        ref_pt = pt
        sol_data = np.array([np.empty_like(res_data[0])])
        while res_data.shape[0]:
            seg, i = self.NearestSeg(res_data, ref_pt)
            sol_data = np.append(sol_data, np.array([seg]), axis = 0)
            res_data = np.delete(res_data, i, axis = 0)
            ref_pt = seg['seg'][1]
        return sol_data[1:], res_data

    def MakePrpDict(self, prp_idx, mtext, lname):
        prp_dict={'feed':200,
                      'ref_coord':np.zeros(3),
                      'power':0,
                      'angle':0,
                      'radius':np.zeros(3),
                      'cut_dir':'cw',
                      'lname':lname,
                      'split':None}

        start_pt = np.zeros(3)

        for text_obj in mtext:
            text = text_obj.get_text()

            d_feed    = re.findall('feed\s*=\s*([\.\d]+)', text)
            d_power   = re.findall('power\s*=\s*([\.\d]+)', text)
            d_angle   = re.findall('angle\s*=\s*([\-\.\d]+)', text)
            d_radius  = re.findall('radius\s*=\s*([\-\.\d]+)', text)
            d_cut_dir = re.findall('cut_dir\s*=\s*(c?cw).*', text)
            d_split   = re.findall('split\s*=\s*([\d]+).*', text)
            d_coord   = re.findall('.*coord_0.*', text, re.IGNORECASE)
            d_start   = re.findall('.*start.*', text, re.IGNORECASE)

            if d_feed:    prp_dict['feed']     = np.float(d_feed[0])
            if d_power:   prp_dict['power']    = np.float(d_power[0])
            if d_angle:   prp_dict['angle']    = np.float(d_angle[0])
            if d_split:   prp_dict['split']    = np.int(d_split[0])
            if d_radius:  prp_dict['radius']   = np.array([0,0,np.float(d_radius[0])])
            if d_cut_dir: prp_dict['cut_dir']  = d_cut_dir
            if d_coord:
                prp_dict['ref_coord']= np.array(text_obj.get_dxf_attrib('insert'))

            if d_start:
                start_pt= np.array(text_obj.get_dxf_attrib('insert'))
                self.prp_dict['glob'].update({'start': start_pt})

            self.prp_dict['loc'].update({prp_idx: prp_dict})

    def AddSct( self, sname):
        print('section name', sname)
        self.model[sname] = {0:{}, 1:{}}
        print(self.model)

    def AddSeg( self, lname):
        new_seg = np.array([], dtype = self.seg_dt)
        prp_idx = hash(lname)
        modelspace = self.dwg.modelspace()
        lines = modelspace.query('LINE[layer=="{}"]'.format(lname))
        arcs = modelspace.query('ARC[layer=="{}"]'.format(lname))
        ltext = modelspace.query('MTEXT[layer=="{}"]'.format(lname))

        if lines:
            seg_lin = np.stack([np.array([ ( np.round(( line.dxf.start, line.dxf.end ),4),prp_idx)], dtype = self.seg_dt) for line in lines], axis = 0)
            new_seg = np.append(new_seg, seg_lin)

        if arcs:
            n = 10;
            min_len = 0.001
            tol = 5

            for var in arcs:
                O =      var.dxf.center
                R =      var.dxf.radius
                angl_1 = var.dxf.start_angle * np.pi / 180
                angl_2 = var.dxf.end_angle * np.pi / 180

                if angl_2 >= angl_1:
                    angl_list = np.linspace(angl_1, angl_2, n)
                else:
                    angl_list = np.linspace(angl_1, angl_2 + 2 * np.pi, n)

                arc_len = R * np.absolute(angl_2 - angl_1)

                if arc_len / n < min_len:
                    n = np.max(np.int(arc_len / min_len), 3)

                ARC_knots_list = []
                for i, angl in enumerate(angl_list):
                    # print(i)
                    ARC_knots_list.append((O[0] + R * np.cos(angl), O[1] + R * np.sin(angl), O[2]))

                for i in range(n - 1):
                    # print(ARC_knots_list[i:i + 2])
                    seg_lin = np.stack([np.array([ ( np.round(( ARC_knots_list[i:i + 2][0], ARC_knots_list[i:i + 2][1] ),4),prp_idx)], dtype = self.seg_dt) for line in lines], axis = 0)
                    # print(( np.round(( ARC_knots_list[i:i + 2][0], ARC_knots_list[i:i + 2][1] ),4)))
                    new_seg = np.append(new_seg, seg_lin)

        self.MakePrpDict(prp_idx, ltext, lname)
        new_seg = self.MakeSplit(prp_idx, new_seg)
        self.seg_arr = np.append(self.seg_arr, new_seg)



    def Seg2Prof(self):
        print('')

    def ApplyTransformations(self):
#apply coord transform
        main_key = list(self.prp_dict['loc'].keys())[0]
        coord_0 = self.prp_dict['loc'][main_key]['ref_coord']
        radius_0 = self.prp_dict['loc'][main_key]['radius']
        self.seg_arr['seg'] -= coord_0
        self.seg_arr['seg'] += radius_0
        self.prp_dict['glob']['start'] -=coord_0

    def MakeChain(self, seg_arr = None):

        if seg_arr is None:
            pt = self.prp_dict['glob']['start']
            buf_data = np.unique(self.OrderSegArr(self.seg_arr), axis = 0)

        else:
            pt = np.zeros(3)
            buf_data = np.unique(self.OrderSegArr(seg_arr), axis = 0)

        chain, rest = self.SortSegArr(buf_data, pt)
        self.seg_sorted = chain
        return chain

    def MakeSplit(self, prp_id, seg_arr):

        if 'split' in self.prp_dict['loc'][prp_id].keys():
            split = self.prp_dict['loc'][prp_id]['split']

            if split:
                chain_arr = self.MakeChain(seg_arr)
                v = chain_arr['seg'][:,1] - chain_arr['seg'][:,0]

                l_norm = np.linalg.norm(v, axis=1)
                l_cs = np.hstack((0, np.cumsum(l_norm)))

                n_seg = np.linspace(0, np.sum(l_norm), split + 1)
                x=np.interp(n_seg, l_cs, np.hstack((chain_arr['seg'][:,0,0], chain_arr['seg'][-1,1,0])))
                y=np.interp(n_seg, l_cs, np.hstack((chain_arr['seg'][:,0,1], chain_arr['seg'][-1,1,1])))
                z=np.interp(n_seg, l_cs, np.hstack((chain_arr['seg'][:,0,2], chain_arr['seg'][-1,1,2])))

                arr_buff0 = np.column_stack((x, y, z))

                arr_buff1= np.roll(arr_buff0, -1, axis=0)
                arr_buff  = np.stack((arr_buff0, arr_buff1),axis=1)[:-1]
                prop_buff = np.ones(arr_buff.shape[0], dtype=np.int) * prp_id

                return np.stack([np.array( ( seg, prp), dtype = self.seg_dt) for seg, prp in zip(arr_buff, prop_buff)], axis = 0)


        return seg_arr




    def PlotChain(self, mode = '2D'):
        if '2D' in mode:
            poly_arr = self.SegArr2Poly(self.seg_sorted['seg'])

            x_arr = poly_arr[:,0]
            y_arr = poly_arr[:,1]

            fig = plt.figure()
            ax = fig.gca()
            ax.plot(x_arr,y_arr)
        plt.ion()
        plt.show()


class AxisProfile():
    def SegArr2Poly(self, seg):
        return np.vstack([seg[:,0], seg[-1,1]])

    def __init__(self):
        self.axis_profile={0:{}}

    def Add2Axis(self, axis_num, chain):
        self.axis_profile[axis_num] = chain

    def Plot(self, mode = '2D'):
        if '2D' in mode:
            fig = plt.figure()
            ax = fig.gca()
            for axis_num ,s in self.axis_profile.items():
                poly_arr = self.SegArr2Poly(s.seg_sorted['seg'])
                x_arr = poly_arr[:,0]
                y_arr = poly_arr[:,1]
                ax.plot(x_arr,y_arr, label = '{}'.format(axis_num))
            plt.legend()
            plt.show()

        if '3D' in mode:
            fig = plt.figure()
            ax = fig.gca(projection='3d')

            for axis_num ,s in self.axis_profile.items():
                poly_arr = self.SegArr2Poly(s.seg_sorted['seg'])
                x_arr = poly_arr[:,0]
                y_arr = poly_arr[:,1]
                z_arr = poly_arr[:,2]
                ax.plot(x_arr, y_arr, z_arr)

            plt.show()

class ModelProfile():
    def transform_pt(self, p, ref_coord = np.zeros(3), r=np.zeros(3), angle=0):
        # def M(axis, theta):
        #     return expm(np.cross(np.eye(3), axis/norm(axis)*np.radians(angle)))
        def M(axis, theta):
            """
            Return the rotation matrix associated with counterclockwise rotation about
            the given axis by theta radians.
            """
            axis = np.asarray(axis)
            axis = axis/np.sqrt(np.dot(axis, axis))
            a = np.cos(theta/2.0)
            b, c, d = -axis*np.sin(theta/2.0)
            aa, bb, cc, dd = a*a, b*b, c*c, d*d
            bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
            return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                             [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                             [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

        axis = np.array([0,0.1,0])
        M0 = M(axis, np.radians(angle))
        p_abs = p+r-ref_coord
        return np.dot(M0, p_abs) #+ np.dot(M0, r)

    def transform_chain(self, chain, prp_arr, prp_dict):
        chain_tr = []
        # print(chain)
        for pt, idx in zip(chain, prp_arr):
            chain_tr.append(self.transform_pt(pt, angle=prp_dict[idx]['angle']))
            # chain_tr.append(self.transform_pt(pt,ref_coord = prp_dict[idx]['ref_coord'], angle=prp_dict[idx]['angle']))
            # chain_tr.append(self.transform_pt(pt,ref_coord = prp_dict[idx]['ref_coord'], angle=prp_dict[idx]['angle'], r=prp_dict[idx]['radius']))
        return np.array(chain_tr)


    def SegArr2Poly(self, seg):
        return np.vstack([seg[:,0], seg[-1,1]])

    def __init__(self):
        self.model_profile={0:{}}

    def Add2Prof(self, prof_num, profile):
        self.model_profile[prof_num] = profile

    def Plot(self, mode = '3D'):

        if '3D' in mode:
            fig = plt.figure()
            ax = fig.gca(projection='3d')

            for prof_num, ap in self.model_profile.items():
                for sec_num, s in ap.axis_profile.items():

                    chain_arr = self.SegArr2Poly(s.seg_sorted['seg'])
                    prp_arr = np.hstack([s.seg_sorted['prp'],s.seg_sorted['prp'][-1]])
                    poly_arr = self.transform_chain(chain_arr, prp_arr, s.prp_dict['loc'])

                    x_arr = poly_arr[:,0]
                    y_arr = poly_arr[:,1]
                    z_arr = poly_arr[:,2]

                    ax.plot(x_arr, y_arr, z_arr)

            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            plt.show()

class CuttingSpace():
    def transform_pt(self, p, ref_coord = np.zeros(3), r=np.zeros(3), angle=0):
        # def M(axis, theta):
        #     return expm(np.cross(np.eye(3), axis/norm(axis)*np.radians(angle)))
        def M(axis, theta):
            """
            Return the rotation matrix associated with counterclockwise rotation about
            the given axis by theta radians.
            """
            axis = np.asarray(axis)
            axis = axis/np.sqrt(np.dot(axis, axis))
            a = np.cos(theta/2.0)
            b, c, d = -axis*np.sin(theta/2.0)
            aa, bb, cc, dd = a*a, b*b, c*c, d*d
            bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
            return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                             [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                             [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

        axis = np.array([0,0.1,0])
        M0 = M(axis, np.radians(angle))
        p_abs = p+r-ref_coord
        return np.dot(M0, p_abs) #+ np.dot(M0, r)

    def transform_chain(self, chain, prp_arr, prp_dict):
        chain_tr = []

        for pt, idx in zip(chain, prp_arr):
            chain_tr.append(self.transform_pt(pt, angle=prp_dict[idx]['angle']))

        return np.array(chain_tr)


    def SegArr2Poly(self, seg):
        return np.vstack([seg[:,0], seg[-1,1]])

    def line_plane_intersect(self, L0, L1, D0, n0):
        """Calculate 3D line-plane intersection point
        Parameters:
            L0, L1 - points on the line
            D0 - point on the plane
            n0 - plane normal
        Returns:
            P - intersetion point. If the line is parallel to the plane,
            function returns [inf, inf, inf]
        """
        l = (L1 - L0)
        n_norm = n0/np.linalg.norm(n0)
        dot_prod = np.dot(l, n_norm)

        if dot_prod == 0:
            P = np.array([np.inf,np.inf,np.inf])
        else:
            d = np.dot((D0 - L0),n_norm)/dot_prod
            P = L0 + d * l
        return P

    def make_projection(self, ch1, ch2, p1, proj_dir):

        proj_ch1 = []

        for var1, var2 in zip(ch1, ch2):
            # print(var1, var2)
            proj_ch1.append(self.line_plane_intersect(var2, var1, np.array([0,0,p1]), np.array([0,0, proj_dir])))

        return np.array(proj_ch1)


    def __init__(self, conf):
        self.conf = conf

    def Add2Cut(self, md):
        self.md = md

    def Plot(self, mode = '3D'):

        rt = self.conf['RTable_loc']
        zs = self.conf['Z_span']
        XYz = -rt
        UVz = -rt + zs

        if '3D' in mode:
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            model_profile = self.md.model_profile
            # print(dir(model_profile))
            for prof_num, ap in model_profile.items():

                if len(ap.axis_profile.keys())==2:

                    print(ap.axis_profile.keys())

                    s = ap.axis_profile
                    # print(s)
                    chain_arr0 = self.SegArr2Poly(s[0].seg_sorted['seg'])
                    chain_arr1 = self.SegArr2Poly(s[1].seg_sorted['seg'])

                    prp_arr0 = np.hstack([s[0].seg_sorted['prp'],s[0].seg_sorted['prp'][-1]])
                    prp_arr1 = np.hstack([s[1]
                    .seg_sorted['prp'],s[1].seg_sorted['prp'][-1]])

                    poly_arr0 = self.transform_chain(chain_arr0, prp_arr0, s[0].prp_dict['loc'])
                    poly_arr1 = self.transform_chain(chain_arr1, prp_arr1, s[1].prp_dict['loc'])

                    poly_arr0proj = self.make_projection(poly_arr0, poly_arr1, XYz, 1)
                    poly_arr1proj = self.make_projection(poly_arr0, poly_arr1, UVz, 1)

                    x_arr0 = poly_arr0[:,0]
                    y_arr0 = poly_arr0[:,1]
                    z_arr0 = poly_arr0[:,2]
                    x_arr1 = poly_arr1[:,0]
                    y_arr1 = poly_arr1[:,1]
                    z_arr1 = poly_arr1[:,2]

                    x_arr0proj = poly_arr0proj[:,0]
                    y_arr0proj = poly_arr0proj[:,1]
                    z_arr0proj = poly_arr0proj[:,2]
                    x_arr1proj = poly_arr1proj[:,0]
                    y_arr1proj = poly_arr1proj[:,1]
                    z_arr1proj = poly_arr1proj[:,2]

                    ax.plot(x_arr0, y_arr0, z_arr0)
                    ax.plot(x_arr1, y_arr1, z_arr1)

                    ax.plot(x_arr0proj, y_arr0proj, z_arr0proj, color='k')
                    ax.plot(x_arr1proj, y_arr1proj, z_arr1proj, color='k')

                else:
                    s = ap.axis_profile[0]
                    chain_arr = self.SegArr2Poly(s.seg_sorted['seg'])
                    prp_arr0 = np.hstack([s.seg_sorted['prp'],s.seg_sorted['prp'][-1]])
                    poly_arr = self.transform_chain(chain_arr, prp_arr0, s.prp_dict['loc'])

                    x_arr = poly_arr[:,0]
                    y_arr = poly_arr[:,1]
                    z_arr = poly_arr[:,2]

                    ax.plot(x_arr, y_arr, z_arr)

            ax.plot(np.array([0,100]), np.zeros(2), np.full((2),XYz))
            ax.plot(np.array([0,100]), np.zeros(2), np.full((2),UVz))
            ax.plot(np.zeros(2), np.array([0,100]), np.full((2),XYz))
            ax.plot(np.zeros(2), np.array([0,100]), np.full((2),UVz))

            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            plt.show()


    def cut_projection2gcode(self, arr0proj, arr1proj, prp_arr0, prp_arr1, prp_dict0, prp_dict1, gcode_conf, machine_conf):
        gcode=[]
        speed_old = 0
        heat_old = 0
        lname_old=''

        for i, (cp1, cp2, c_prop1) in enumerate(zip(arr0proj, arr1proj, prp_arr0)):

            lname_new=prp_dict0['loc'][c_prop1]['lname']
            speed_new=prp_dict0['loc'][c_prop1]['feed']
            heat_new= prp_dict0['loc'][c_prop1]['power']
            angle =   prp_dict0['loc'][c_prop1]['angle']
            # ref_coord = prp_dict0[c_prop1]['ref_coord']
            if lname_new != lname_old:
                line = '({1})'.format(gcode_conf['spindle'],
                            lname_new)
                gcode.append(line)
                lname_old = lname_new

            if heat_new != heat_old:
                line = '{0[0]} {0[1]}{1}'.format(gcode_conf['spindle'],
                                        heat_new)
                gcode.append(line)
                heat_old = heat_new

            if speed_new != speed_old:
                line = '{0} {1}'.format(gcode_conf['speed'],
                                        speed_new)
                gcode.append(line)
                speed_old = speed_new

            if i==0:
                line = 'G1 {1[4]} {0[2]:<8.2f} '.format([cp1, cp2, angle], machine_conf['ax'])
                gcode.append(line)

            line = 'G1 {1[0]} {0[0][0]:<8.2f} ' \
                      '{1[1]} {0[0][1]:<8.2f} ' \
                      '{1[2]} {0[1][0]:<8.2f} ' \
                      '{1[3]} {0[1][1]:<8.2f} ' \
                      '{1[4]} {0[2]:<8.2f} '.format([cp1, cp2, angle], machine_conf['ax'])

            gcode.append(line)

        gcode.append('{0[2]}'.format(gcode_conf['spindle']))
        return gcode

    def proj_stats(self, master_io_path, master_io_prop, master_prop_dict):
        '''
        TODO: add rotation time!
        '''
        def ct_sect2(sect_arr):
            p = sect_arr - np.roll(sect_arr,-1, axis=0)
            l_arr = np.linalg.norm( p, axis=1)
            return l_arr

        def ct_speed2(l_arr, prop, prop_dict):
            # print(prop)
            # print('mmmmmmmmmm',prop_dict['loc'][prop[0]]['feed'])
            if l_arr.size:
                sect_speed = ct_sect2(l_arr) * 60 / np.array([prop_dict['loc'][key]['feed'] for key in prop])
            else:
                sect_speed = 0
            return np.sum(sect_speed)

        def ct_len_2(l_arr):
            # print(l_arr)
            if l_arr.size:
                res = np.sum(ct_sect2(l_arr))
            else:
                res = 0
            return res

        ct_len_list = []
        ct_time_list = []
        # print('{:-^79}'.format('CUTTING STATS'))
        # print(master_io_path)
        io_ct_len = ct_len_2(master_io_path)
        # print('ct len',io_ct_len)
        # io_ct_len = 0
        # io_ct_speed =0
        io_ct_speed = ct_speed2(master_io_path, master_io_prop, master_prop_dict)
        # print('layer:  {}'.format(master_prop_dict[0]['layer']))
        return io_ct_len, io_ct_speed


    def SaveGcode(self, gcode_fname):
        gcode=[]
        rt = self.conf['RTable_loc']
        zs = self.conf['Z_span']
        XYz = -rt
        UVz = -rt + zs
        machine_conf={
    #axA master column
        'ax':['X','Y', 'U', 'V', 'B'],
    #distance between columns.
        'AB_dist':0.45,
    # distance between rotary table and column A
        'AT_dist':0.225}
        gcode_conf={'comment':['(',')'],
                    'spindle':['M3','S','M5'],
                    'speed':'F',
                    'end':'M2'}

        ct_len_list = []
        ct_time_list = []

        tot_ct_l = 0
        tot_ct_t = 0

        model_profile = self.md.model_profile

        # print(ss)
        print('{:-^79}'.format('GCODE'))

        # for i, var1 in enumerate(ss):

        for prof_num, ap in model_profile.items():

            s = ap.axis_profile

            if len(ap.axis_profile.keys())==2:

                # print(ap.axis_profile.keys())

                # print(s)
                chain_arr0 = self.SegArr2Poly(s[0].seg_sorted['seg'])
                chain_arr1 = self.SegArr2Poly(s[1].seg_sorted['seg'])

                prp_arr0 = np.hstack([s[0].seg_sorted['prp'],s[0].seg_sorted['prp'][-1]])
                prp_arr1 = np.hstack([s[1].seg_sorted['prp'],s[1].seg_sorted['prp'][-1]])

                poly_arr0 = chain_arr0
                poly_arr1 = chain_arr1

                poly_arr0proj = self.make_projection(poly_arr0, poly_arr1, XYz, 1)
                poly_arr1proj = self.make_projection(poly_arr0, poly_arr1, UVz, 1)

                arr0proj = poly_arr0proj[:,:2]
                arr1proj = poly_arr1proj[:,:2]

                prp_dict0 = s[0].prp_dict
                prp_dict1 = s[1].prp_dict
                # print(s[0].prp_dict)
                # print(s[1].prp_dict)

            else:
                chain_arr0 = self.SegArr2Poly(s[0].seg_sorted['seg'])

                prp_arr0 = np.hstack([s[0].seg_sorted['prp'],s[0].seg_sorted['prp'][-1]])
                prp_arr1 = np.hstack([s[0].seg_sorted['prp'],s[0].seg_sorted['prp'][-1]])

                poly_arr0 = chain_arr0

                arr0proj=poly_arr0[:,:2]
                arr1proj=poly_arr0[:,:2]

                prp_dict0 = s[0].prp_dict
                prp_dict1 = s[0].prp_dict
                # print(s[0].prp_dict)

        #
            # gcode.append('{0[0]}layers: {1[0]} {1[1]}{0[1]}'.format(gcode_conf['comment'],[prop_dict1[0]['layer'], prop_dict2[0]['layer']]))
            # gcode.append('{0[0]}sequence :{1}{0[1]}'.format(gcode_conf['comment'],i))
        #
        #
            ct_l, ct_t = self.proj_stats(arr0proj, prp_arr0, prp_dict0)
        # #
        #     print('{:-^79}'.format(prop_dict1[0]['layer']))
        #     print('{0:8s}{1:10.0f}mm'.format('length:', ct_l))
        #     print('{0:8s}{1:10.0f}s'.format('time:', ct_t))
        # #
            tot_ct_l += ct_l
            tot_ct_t += ct_t
        #
        #
        #

            gcode += self.cut_projection2gcode(arr0proj, arr1proj, prp_arr0, prp_arr1, prp_dict0, prp_dict1, gcode_conf, machine_conf)

        print('{:-^79}'.format('TOTAL'))
        print('{0:8s}{1:10.0f}mm'.format('length:', tot_ct_l))
        print('{0:8s}{1:10.0f}s'.format('time:', tot_ct_t))

        gcode.append('{0}'.format(gcode_conf['end']))
        gcode.append('%')



        with open(gcode_fname,'w') as f:
            for line in gcode:
                f.write(line+'\n')
