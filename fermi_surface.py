import re
import numpy as np
import plotly.graph_objects as go
from pymatgen.io import vasp

class fermi_surface():
    def __init__(self, filename, vasprun='vasprun.xml', fermi=0):
        self.readfile(filename)
        self.fermi=fermi
        self.vasprun = vasp.Vasprun(vasprun)
        
    def readfile(self, filename):
        pattern = re.compile(r"[+-]?(?:[0-9]*[.])?[0-9]+")
        f = open(filename, 'r')
        data = f.readline()
        self.ngrid = np.array(pattern.findall(data), dtype=int)

        data = f.readline()
        self.nmode = int(pattern.findall(data)[0])

        data = f.readline()
        self.nband = int(pattern.findall(data)[0])

        self.rec_latt = np.zeros((3,3))
        for i in range(3):
            data = f.readline()
            self.rec_latt[i]=np.array(pattern.findall(data), dtype=np.float)

        self.energy = np.zeros(np.append(self.nband, self.ngrid))
        for iband in range(self.nband):
            for i in range(self.ngrid[0]):
                for j in range(self.ngrid[1]):
                    for k in range(self.ngrid[2]):
                        data = f.readline()
                        self.energy[iband, i,j,k] = float(pattern.findall(data)[0])

        self.color = np.zeros(np.append(self.nband, self.ngrid))
        for iband in range(self.nband):
            for i in range(self.ngrid[0]):
                for j in range(self.ngrid[1]):
                    for k in range(self.ngrid[2]):
                        data = f.readline()
                        self.color[iband, i,j,k] = float(pattern.findall(data)[0])

    def interp(self, e1, e2, c1, c2):
        if np.abs(e1-e2)<1e-12:
            return (c1+c2)/2
        return e1/(e1-e2)*c2+e2/(e2-e1)*c1

    def add_triangle(self, kvecs, eigc):
        kvecs = np.array(kvecs)
        self.tx.append(kvecs.T[0])
        self.ty.append(kvecs.T[1])
        self.tz.append(kvecs.T[2])
        self.tc.append(eigc)

    def patch_triangle(self):
        corner = [[ 0, 1, 3, 7 ],
          [ 0, 1, 5, 7 ],
          [ 0, 2, 3, 7 ],
          [ 0, 2, 6, 7 ],
          [ 0, 4, 5, 7 ],
          [ 0, 4, 6, 7 ]]

        dr = [[0, 0, 0],
            [0, 0, 1],
            [0, 1, 0],
            [0, 1, 1],
            [1, 0, 0],
            [1, 0, 1],
            [1, 1, 0],
            [1, 1, 1]]

        self.tx = []
        self.ty = []
        self.tz = []
        self.tc = []

        for ib in range(self.nband):
            for ix in range(self.ngrid[0]):
                for iy in range(self.ngrid[1]):
                    for iz in range(self.ngrid[2]):
                        eig = np.zeros(8)
                        eig_c = np.zeros(8)
                        kvec = np.zeros((8,3))
                        ngrid = self.ngrid
                        for i in range(8):
                            eig[i] = self.energy[ib][(ix+dr[i][0])%ngrid[0]][(iy+dr[i][1])%ngrid[1]][(iz+dr[i][2])%ngrid[2]] - self.fermi
                            eig_c[i] = self.color[ib][(ix+dr[i][0])%ngrid[0]][(iy+dr[i][1])%ngrid[1]][(iz+dr[i][2])%ngrid[2]]
                            kvec[i] = (np.array([ix+dr[i][0], iy+dr[i][1], iz+dr[i][2]]) / ngrid).dot(self.rec_latt)
                        for it in range(6):
                            eig_t = eig[corner[it]]
                            eig_sort = np.sort(eig_t)
                            order = np.argsort(eig_t)
                            kvec_sort = kvec[corner[it]][order]
                            eig_c_sort = eig_c[corner[it]][order]
                            if eig_sort[0]>0:
                                pass
                            elif eig_sort[1]>0:
                                k1 = self.interp(eig_sort[0], eig_sort[1], kvec_sort[0], kvec_sort[1])
                                k2 = self.interp(eig_sort[0], eig_sort[2], kvec_sort[0], kvec_sort[2])
                                k3 = self.interp(eig_sort[0], eig_sort[3], kvec_sort[0], kvec_sort[3])
                                c1 = self.interp(eig_sort[0], eig_sort[1], eig_c_sort[0], eig_c_sort[1])
                                c2 = self.interp(eig_sort[0], eig_sort[2], eig_c_sort[0], eig_c_sort[2])
                                c3 = self.interp(eig_sort[0], eig_sort[3], eig_c_sort[0], eig_c_sort[3])
                                self.add_triangle([k1, k2, k3], [c1, c2, c3])
                            elif eig_sort[2]>0:
                                k2 = self.interp(eig_sort[0], eig_sort[3], kvec_sort[0], kvec_sort[3])
                                k3 = self.interp(eig_sort[1], eig_sort[3], kvec_sort[1], kvec_sort[3])
                                c1 = self.interp(eig_sort[0], eig_sort[2], eig_c_sort[0], eig_c_sort[2])
                                c2 = self.interp(eig_sort[0], eig_sort[3], eig_c_sort[0], eig_c_sort[3])
                                k1 = self.interp(eig_sort[0], eig_sort[2], kvec_sort[0], kvec_sort[2])
                                c3 = self.interp(eig_sort[1], eig_sort[3], eig_c_sort[1], eig_c_sort[3])
                                self.add_triangle([k1, k2, k3], [c1, c2, c3])

                                k1 = self.interp(eig_sort[0], eig_sort[2], kvec_sort[0], kvec_sort[2])
                                k2 = self.interp(eig_sort[1], eig_sort[2], kvec_sort[1], kvec_sort[2])
                                k3 = self.interp(eig_sort[1], eig_sort[3], kvec_sort[1], kvec_sort[3])
                                c1 = self.interp(eig_sort[0], eig_sort[2], eig_c_sort[0], eig_c_sort[2])
                                c2 = self.interp(eig_sort[1], eig_sort[2], eig_c_sort[1], eig_c_sort[2])
                                c3 = self.interp(eig_sort[1], eig_sort[3], eig_c_sort[1], eig_c_sort[3])
                                self.add_triangle([k1, k2, k3], [c1, c2, c3])
                            elif eig_sort[3]>0:
                                k1 = self.interp(eig_sort[0], eig_sort[3], kvec_sort[0], kvec_sort[3])
                                k2 = self.interp(eig_sort[1], eig_sort[3], kvec_sort[1], kvec_sort[3])
                                k3 = self.interp(eig_sort[2], eig_sort[3], kvec_sort[2], kvec_sort[3])
                                c1 = self.interp(eig_sort[0], eig_sort[3], eig_c_sort[0], eig_c_sort[3])
                                c2 = self.interp(eig_sort[1], eig_sort[3], eig_c_sort[1], eig_c_sort[3])
                                c3 = self.interp(eig_sort[2], eig_sort[3], eig_c_sort[2], eig_c_sort[3])
                                self.add_triangle([k1, k2, k3], [c1, c2, c3])
    
    def get_boundary(self):
        pass
    
    def draw_boundary(self):
        pass

    def move_triangle(self, kvecs, eigc):
        center = []
        shift = []
        for ikx in np.arange(-1, 2):
            for iky in np.arange(-1, 2):
                for ikz in np.arange(-1, 2):
                    if not all([ikx==0, iky==0, ikz==0]):
                        shift.append([ikx, iky, ikz])
                        center.append(np.array([ikx, iky, ikz]).dot(self.rec_latt))
        
        thre = 1e-8
        kvecs = np.array(kvecs)
        if np.linalg.norm(kvecs[0]-kvecs[1])<thre or np.linalg.norm(kvecs[0]-kvecs[2])<thre or np.linalg.norm(kvecs[1]-kvecs[2])<thre:
            return
        eigc = np.array(eigc)
        
        for i in range(len(center)):
            dist = np.linalg.norm(center[i])
            prod = np.dot(kvecs, center[i])/dist
            prod_sort = np.sort(prod)
            order = np.argsort(prod)
            kvecs_sort = kvecs[order]
            eigc_sort = eigc[order]
            a = [[(0.5*dist-prod_sort[j])/(prod_sort[i]-prod_sort[j]) for j in range(3)] for i in range(3)] 
            if prod_sort[0]>dist/2-thre:
                self.move_triangle(kvecs-center[i], eigc)
                return
            elif prod_sort[1]>dist/2+thre:
                x0 = kvecs_sort[0]
                c0 = eigc_sort[0]
                x1 = kvecs_sort[0] * a[0][1] + kvecs_sort[1] * a[1][0]
                c1 = eigc_sort[0] * a[0][1] + eigc_sort[1] * a[1][0]
                x2 = kvecs_sort[0] * a[0][2] + kvecs_sort[2] * a[2][0]
                c2 = eigc_sort[0] * a[0][2] + eigc_sort[2] * a[2][0]
                self.move_triangle([x0, x1, x2], [c0, c1, c2])
                self.move_triangle([kvecs_sort[1], kvecs_sort[2], x1], [eigc_sort[1], eigc_sort[2], c1])
                self.move_triangle([kvecs_sort[2], x1, x2], [eigc_sort[2], c1, c2])
                return
            elif prod_sort[2]>dist/2+thre:
                x0 = kvecs_sort[0]
                c0 = eigc_sort[0]
                x1 = kvecs_sort[1]
                c1 = eigc_sort[1]
                x2 = kvecs_sort[0] * a[0][2] + kvecs_sort[2] * a[2][0]
                c2 = eigc_sort[0] * a[0][2] + eigc_sort[2] * a[2][0]
                x3 = kvecs_sort[1] * a[1][2] + kvecs_sort[2] * a[2][1]
                c3 = eigc_sort[1] * a[1][2] + eigc_sort[2] * a[2][1]
                self.move_triangle([x0, x1, x2], [c0, c1, c2])
                self.move_triangle([kvecs_sort[1], x2, x3],[eigc_sort[1], c2, c3])
                self.move_triangle([kvecs_sort[2], x2+thre, x3+thre],[eigc_sort[2], c2, c3])
                return
        self.triangles.append(kvecs)
        self.triangles_c.append(eigc)

    def move_triangle_bz(self):
        self.triangles = []
        self.triangles_c = []
        for s in range(len(self.tx)):
            self.move_triangle([[self.tx[s][0], self.ty[s][0], self.tz[s][0]],
                            [self.tx[s][1], self.ty[s][1], self.tz[s][1]],
                            [self.tx[s][2], self.ty[s][2], self.tz[s][2]]], 
                            [self.tc[s][0], self.tc[s][1], self.tc[s][2]])

    def draw(self, fermi=8.34):
        self.fermi=fermi
        trace = []
        BZ = self.vasprun.final_structure.lattice.get_brillouin_zone()
        for boundary in BZ:
            x, y, z = np.append(boundary, [boundary[0]], axis=0).T
            trace.append(go.Scatter3d(
                x=x, y=y, z=z, mode='lines', legendgroup='BZ', name='1st Brillouin Zone',
                line=dict(color='darkblue', width=2),
                showlegend=False,
                hoverinfo='none'
            ))
        trace[0]['showlegend']=True
        trace.append(go.Mesh3d(
                x=np.array(self.triangles).reshape(-1, 3)[:,0],
                y=np.array(self.triangles).reshape(-1, 3)[:,1],
                z=np.array(self.triangles).reshape(-1, 3)[:,2],
                colorbar_title='$log <b^2>$',
                colorscale=[[0, 'gold'],
                            [0.5, 'mediumturquoise'],
                            [1, 'magenta']],
                # Intensity of each vertex, which will be interpolated and color-coded
                intensity=np.array(self.triangles_c).reshape(-1),
                # i, j and k give the vertices of triangles
                # here we represent the 4 triangles of the tetrahedron surface
                i=np.arange(0,len(self.triangles)*3,3),
                j=np.arange(1,len(self.triangles)*3,3),
                k=np.arange(2,len(self.triangles)*3,3),
                name='y',
                showscale=True))

        fig = go.Figure(data=trace)
        fig.update_layout(
            width=800,
            height=700,
            autosize=False,
            scene=dict(
                camera=dict(
                    up=dict(
                        x=0,
                        y=0,
                        z=1
                    ),
                    eye=dict(
                        x=1,
                        y=1,
                        z=0.8,
                    )
                ),
                aspectratio = dict( x=1, y=1, z=1 ),
                aspectmode = 'manual',
                camera_projection=dict(type='perspective'),
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False))
        )
        fig.update_layout(legend_orientation="h")
        #fig.show()
        return fig

    def draw_hist(self):
        fig = go.Figure(data = [go.Histogram(x=np.array(self.triangles_c).reshape(-1), nbinsx=30)])
        fig.update_layout(width=800, height=700)
        return fig