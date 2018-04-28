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
        print('entered sort seg arr')
        res_data = buf_data
        ref_pt = pt
        sol_data = np.array([np.empty_like(res_data[0])])
        while res_data.shape[0]:
            seg, i = self.NearestSeg(res_data, ref_pt)
            sol_data = np.append(sol_data, np.array([seg]), axis = 0)
            res_data = np.delete(res_data, i, axis = 0)
            ref_pt = seg['seg'][1]
        return sol_data[1:], res_data

    def MakePrpDict(self, prp_idx,  mtext):
        prp_dict={'feed':200,
                      'ref_coord':np.zeros(3),
                      'power':0,
                      'angle':0,
                      'radius':np.zeros(3),
                      'cut_dir':'cw',
                      # 'layer':layer,
                      'split':None}

        start_pt = np.zeros(3)

        for text_obj in mtext:
            text = text_obj.get_text()

            d_feed    = re.findall('feed\s*=\s*([\.\d]+)', text)
            d_power   = re.findall('power\s*=\s*([\.\d]+)', text)
            d_angle   = re.findall('angle\s*=\s*([\-\.\d]+)', text)
            d_radius  = re.findall('radius\s*=\s*([\-\.\d]+)', text)
            d_cut_dir = re.findall('cut_dir.*=.*(c?cw).*', text)
            d_split   = re.findall('split.*=.*([\d]+).*', text)
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
        prp_idx = hash(lname)

        modelspace = self.dwg.modelspace()
        lines = modelspace.query('LINE[layer=="{}"]'.format(lname))
        arcs = modelspace.query('ARC[layer=="{}"]'.format(lname))
        ltext = modelspace.query('MTEXT[layer=="{}"]'.format(lname))

        if lines:
            seg_lin = np.stack([np.array([ ( np.round(( line.dxf.start, line.dxf.end ),4),prp_idx)], dtype = self.seg_dt) for line in lines], axis = 0)
            self.seg_arr = np.append(self.seg_arr, seg_lin)

        if arcs:
            n = 10;
            min_len = 0.001
            tol = 5

            for var in arcs:
                # print('dasdasdasdasadsasssssssssssssssssssssssssss')
                # print(var.dxf.start_angle)
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
                    self.seg_arr = np.append(self.seg_arr, seg_lin)

        # self.prp_dict.update(self.MakePrpDict(prp_idx,ltext))
        self.MakePrpDict(prp_idx,ltext)
        # print('------------------------------')
        # print('self prp dict  ',self.prp_dict)

    def PrintList(self):
        # pp.pprint( np.vstack(self.seg_arr) )
        # pp.pprint( self.prp_dict )
        print('')
        # print(dxf_read_2(self.dwg, lname))

    def ApplyTransformations(self):
        #apply coord transform
        main_key = list(self.prp_dict['loc'].keys())[0]
        coord_0 = self.prp_dict['loc'][main_key]['ref_coord']
        radius_0 = self.prp_dict['loc'][main_key]['radius']
        # print(coord_0)
        # print(self.prp_dict['glob']['start'])
        self.seg_arr['seg'] -= coord_0
        self.seg_arr['seg'] += radius_0
        self.prp_dict['glob']['start'] -=coord_0

    def MakeChain(self):
        pt = self.prp_dict['glob']['start']
        buf_data = np.unique(self.OrderSegArr(self.seg_arr), axis = 0)
        r1, r2 = self.SortSegArr(buf_data, pt)
        self.seg_sorted = r1
        return r1
#
    def PlotChain(self, mode = '2D'):
        if '2D' in mode:
            poly_arr = self.SegArr2Poly(self.seg_sorted['seg'])

            x_arr = poly_arr[:,0]
            y_arr = poly_arr[:,1]

            fig = plt.figure()
            ax = fig.gca()
            ax.plot(x_arr,y_arr)
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
        print(chain)
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
        print(chain)
        for pt, idx in zip(chain, prp_arr):
            chain_tr.append(self.transform_pt(pt, angle=prp_dict[idx]['angle']))
            # chain_tr.append(self.transform_pt(pt,ref_coord = prp_dict[idx]['ref_coord'], angle=prp_dict[idx]['angle']))
            # chain_tr.append(self.transform_pt(pt,ref_coord = prp_dict[idx]['ref_coord'], angle=prp_dict[idx]['angle'], r=prp_dict[idx]['radius']))
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

    def make_projection(self, chain1, chain2, p1):

        proj_ch1 = []

        for var1, var2 in zip(ch1, ch2):
            proj_ch1.append(self.line_plane_intersect(var2, var1, np.array([0,0,-310]), np.array([0,0,1])))

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

                    s = ap.axis_profile
                    chain_arr0 = self.SegArr2Poly(s.seg_sorted['seg'])
                    chain_arr1 = self.SegArr2Poly(s.seg_sorted['seg'])
                    prp_arr0 = np.hstack([s.seg_sorted['prp'],s.seg_sorted['prp'][-1]])
                    prp_arr1 = np.hstack([s.seg_sorted['prp'],s.seg_sorted['prp'][-1]])
                    poly_arr0 = self.transform_chain(chain_arr, prp_arr, s.prp_dict['loc'])
                    poly_arr1 = self.transform_chain(chain_arr, prp_arr, s.prp_dict['loc'])

                    x_arr0 = poly_arr0[:,0]
                    y_arr0 = poly_arr0[:,1]
                    z_arr0 = poly_arr0[:,2]

                    x_arr1 = poly_arr1[:,0]
                    y_arr1 = poly_arr1[:,1]
                    z_arr1 = poly_arr1[:,2]

                    ax.plot(x_arr0, y_arr0, z_arr0)
                    ax.plot(x_arr1, y_arr1, z_arr1)


                else:
                    s = ap.axis_profile[0]
                    chain_arr = self.SegArr2Poly(s.seg_sorted['seg'])
                    prp_arr = np.hstack([s.seg_sorted['prp'],s.seg_sorted['prp'][-1]])
                    poly_arr = self.transform_chain(chain_arr, prp_arr, s.prp_dict['loc'])

                    x_arr = poly_arr[:,0]
                    y_arr = poly_arr[:,1]
                    z_arr = poly_arr[:,2]

                    ax.plot(x_arr, y_arr, z_arr)

                # for sec_num, s in ap.axis_profile.items():


            ax.plot(np.array([0,100]), np.zeros(2), np.full((2),XYz))
            ax.plot(np.array([0,100]), np.zeros(2), np.full((2),UVz))
            ax.plot(np.zeros(2), np.array([0,100]), np.full((2),XYz))
            ax.plot(np.zeros(2), np.array([0,100]), np.full((2),UVz))
            # ax.plot(x_arr, y_arr, z_arr)
            # ax.plot(x_arr, y_arr, z_arr)
            # ax.plot(x_arr, y_arr, z_arr)

            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            plt.show()
