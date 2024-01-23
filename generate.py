#!/usr/bin/env python

# Copyright (c) 2024, Ronald Römer
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

import math
import os
from functools import reduce

from vtkmodules.vtkCommonCore import vtkIdList, vtkPoints
from vtkmodules.vtkCommonDataModel import vtkPolyData, VTK_POLYGON
from vtkmodules.vtkIOGeometry import vtkSTLWriter

from vtkmodules.vtkFiltersCore import vtkPolyDataNormals, vtkCleanPolyData, vtkTriangleFilter
from vtkmodules.vtkFiltersModeling import vtkLinearExtrusionFilter

# only works with vtkbool#3d0d954

from vtkbool.vtkBool import vtkPolyDataBooleanFilter

class Chimney:
    def __init__(self, cfg):
        self.cfg = { 'e': .15, 'f': .625, 'h': .75, 'l': 1, 'div_a': 1, 'div_b': 1, 'phi': 0, 's': 0 }

        self.cfg.update(cfg)

        assert all( k in self.cfg for k in ['a', 'b', 'seqs', 'count', 't'] )

        assert all( c in [2, 3, 4] for c in sum(map(list, sum(self.cfg['seqs'], [])), []) )

        self.counts_a = [ sum(s_a) for s_a, _ in self.cfg['seqs'] ]
        self.counts_b = [ sum(s_b) for _, s_b in self.cfg['seqs'] ]

        assert self.counts_a.count(self.counts_a[0]) == len(self.counts_a)
        assert self.counts_b.count(self.counts_b[0]) == len(self.counts_b)

        assert self.cfg['s'] >= 0 and self.cfg['s'] <= self.cfg['a']

    @staticmethod
    def transform(line, ang, dx, dy):
        tranformed = []

        for p in line:
            x = math.cos(ang)*p[0]-math.sin(ang)*p[1]
            y = math.sin(ang)*p[0]+math.cos(ang)*p[1]

            tranformed.append((x+dx, y+dy))

        return tranformed

    @staticmethod
    def get_line(w, e, seq):
        x = (w+e)/sum(seq)

        xs = [0]

        for c in seq[:-1]:
            xs.append(xs[-1]+c*x)

        line = []

        for s, l in zip(xs, map(lambda c: c*x, seq)):
            line.extend([(s, 0), (s+l-e, 0), (s+l-e, e), (s+l, e)])

        del line[-3:]

        return line

    @staticmethod
    def get_bricks_line(a, seq_a, b, seq_b, e):
        line_a = Chimney.get_line(a, e, seq_a)
        line_b = Chimney.get_line(b, e, seq_b)

        pl = line_a + Chimney.transform(line_b, math.pi/2, a, 0) + Chimney.transform(line_a, math.pi, a, b) + Chimney.transform(line_b, 3*math.pi/2, 0, b)

        return pl

    @staticmethod
    def append_z(line, z=0):
        return [ (*pt, z) for pt in line ]

    @staticmethod
    def extrude_from_line(line, x=0, y=0, z=0):
        cell = vtkIdList()

        pts = vtkPoints()
        pts.SetDataTypeToDouble()

        [ (pts.InsertNextPoint(pt), cell.InsertNextId(i)) for i, pt in enumerate(line) ]

        pd = vtkPolyData()
        pd.Allocate(1)
        pd.SetPoints(pts)
        pd.InsertNextCell(VTK_POLYGON, cell)

        extr = vtkLinearExtrusionFilter()
        extr.SetInputData(pd)
        extr.SetVector(x, y, z)

        normals = vtkPolyDataNormals()
        normals.SetInputConnection(extr.GetOutputPort())
        normals.AutoOrientNormalsOn()

        return normals

    @staticmethod
    def union(f_a, f_b):
        bf = vtkPolyDataBooleanFilter()
        bf.SetInputConnection(0, f_a.GetOutputPort())
        bf.SetInputConnection(1, f_b.GetOutputPort())

        return bf

    @staticmethod
    def difference(f_a, f_b):
        bf = vtkPolyDataBooleanFilter()
        bf.SetInputConnection(0, f_a.GetOutputPort())
        bf.SetInputConnection(1, f_b.GetOutputPort())
        bf.SetOperModeToDifference()

        return bf

    def export(self, name):
        all_seqs = []

        while len(all_seqs) < self.cfg['count']:
            all_seqs.extend(self.cfg['seqs'])

            del all_seqs[self.cfg['count']:]

        g = self.cfg['f']-self.cfg['e']

        bfs = []

        for i, (s_a, s_b) in enumerate(all_seqs, start=1):
            z = -i*self.cfg['f']

            bricks_line = Chimney.get_bricks_line(self.cfg['a'], s_a, self.cfg['b'], s_b, self.cfg['e'])

            bricks = Chimney.extrude_from_line(Chimney.append_z(bricks_line, z), z=g)

            gap_line = [(self.cfg['e'], self.cfg['e']),
                (self.cfg['a']-self.cfg['e'], self.cfg['e']),
                (self.cfg['a']-self.cfg['e'], self.cfg['b']-self.cfg['e']),
                (self.cfg['e'], self.cfg['b']-self.cfg['e'])]

            gap = Chimney.extrude_from_line(Chimney.append_z(gap_line, z+g), z=self.cfg['e'])

            bfs.append(Chimney.union(bricks, gap))

        bf = reduce(lambda a, b: Chimney.union(a, b), bfs)

        chimney = Chimney.union(bf, Chimney.extrude_from_line(Chimney.append_z([(0, 0), (self.cfg['a'], 0), (self.cfg['a'], self.cfg['b']), (0, self.cfg['b'])], 0), z=self.cfg['h']))

        st_a = (self.cfg['a']+self.cfg['e'])/self.counts_a[0]-self.cfg['e']/2
        st_b = (self.cfg['b']+self.cfg['e'])/self.counts_b[0]-self.cfg['e']/2

        xs = [ st_a+i*(self.cfg['a']-2*st_a)/self.cfg['div_a'] for i in range(self.cfg['div_a']+1) ]
        ys = [ st_b+i*(self.cfg['b']-2*st_b)/self.cfg['div_b'] for i in range(self.cfg['div_b']+1) ]

        bfs_ = [chimney]

        for x0, x1 in zip(xs, xs[1:]):
            for y0, y1 in zip(ys, ys[1:]):
                out = Chimney.extrude_from_line(Chimney.append_z([(x0+st_a, y0+st_b), (x1-st_a, y0+st_b), (x1-st_a, y1-st_b), (x0+st_a, y1-st_b)], self.cfg['h']), z=-5)

                bfs_.append(out)

        chimney_ = reduce(lambda a, b: Chimney.difference(a, b), bfs_)

        off = 1

        line = [(self.cfg['s'], -off, -self.cfg['t']),
            (-off, -off, -self.cfg['t']-math.tan(self.cfg['phi'])*(self.cfg['s']+off)),
            (-off, -off, -self.cfg['count']*self.cfg['f']-self.cfg['h']-off),
            (self.cfg['a']+off, -off, -self.cfg['count']*self.cfg['f']-self.cfg['h']-off),
            (self.cfg['a']+off, -off, -self.cfg['t']-math.tan(self.cfg['phi'])*(self.cfg['a']-self.cfg['s']+off))]

        result = Chimney.difference(chimney_, Chimney.extrude_from_line(line, y=self.cfg['b']+2*off))

        base_line = [(self.cfg['s'], 0, -self.cfg['t']),
            (0, 0, -self.cfg['t']-math.tan(self.cfg['phi'])*self.cfg['s']),
            (0, 0, -self.cfg['t']-self.cfg['l']-math.tan(self.cfg['phi'])*self.cfg['s']),
            (self.cfg['s'], 0, -self.cfg['t']-self.cfg['l']),
            (self.cfg['a'], 0, -self.cfg['t']-self.cfg['l']-math.tan(self.cfg['phi'])*(self.cfg['a']-self.cfg['s'])),
            (self.cfg['a'], 0, -self.cfg['t']-math.tan(self.cfg['phi'])*(self.cfg['a']-self.cfg['s']))]

        if self.cfg['s'] == 0 or self.cfg['s'] == self.cfg['a']:
            del base_line[0]

        base = Chimney.extrude_from_line(base_line, y=self.cfg['b'])

        result_ = Chimney.union(result, base)

        clean = vtkCleanPolyData()
        clean.SetInputConnection(result_.GetOutputPort())

        tf = vtkTriangleFilter()
        tf.SetInputConnection(clean.GetOutputPort())

        writer = vtkSTLWriter()
        writer.SetInputConnection(tf.GetOutputPort())
        writer.SetFileName(name)
        writer.Update()

if __name__ == '__main__':
    cfgs = [
        { 'a': 4.3333, 'b': 6.5, 'seqs': [
            [(2, 2, 2, 2), (3, 4, 2, 3)],
            [(3, 2, 3), (2, 2, 2, 2, 2, 2)],
            [(2, 2, 2, 2), (3, 2, 4, 3)],
            [(3, 2, 3), (2, 2, 2, 2, 2, 2)]
        ], 'div_a': 1, 'div_b': 2, 'count': 18, 's': 4.3333, 't': 6.916667, 'phi': math.pi/4 },
        { 'a': 3.25, 'b': 8.6666, 'seqs': [
            [(3, 3), (2, 2, 2, 2, 2, 2, 2, 2)],
            [(2, 2, 2), (3, 4, 4, 2, 3)],
            [(3, 3), (2, 2, 2, 2, 2, 2, 2, 2)],
            [(2, 2, 2), (3, 2, 4, 4, 3)],
        ], 'div_a': 1, 'div_b': 3, 'count': 13, 's': 3.25, 't': 4.875, 'phi': math.pi/4 },
        { 'a': 10.8333, 'b': 4.3333, 'seqs': [
            [(3, 4, 4, 4, 2, 3), (2, 2, 2, 2)],
            [(2, 2, 2, 2, 2, 2, 2, 2, 2, 2), (3, 2, 3)],
            [(3, 2, 4, 4, 4, 3), (2, 2, 2, 2)],
            [(2, 2, 2, 2, 2, 2, 2, 2, 2, 2), (3, 2, 3)]
        ], 'div_a': 3, 'div_b': 1, 'count': 15, 's': 7.583333, 't': 1.75, 'phi': math.pi/4 }
    ]

    os.makedirs('out', exist_ok=True)

    for i, cfg in enumerate(cfgs):
        Chimney(cfg).export(f'out/chimney{i}.stl')
