#!/usr/bin/env python

import math
from functools import reduce

from vtkmodules.vtkCommonCore import vtkIdList, vtkPoints
from vtkmodules.vtkCommonDataModel import vtkPolyData, VTK_POLYGON
from vtkmodules.vtkIOLegacy import vtkPolyDataWriter
from vtkmodules.vtkIOGeometry import vtkSTLWriter

from vtkmodules.vtkFiltersCore import vtkPolyDataNormals, vtkCleanPolyData
from vtkmodules.vtkFiltersModeling import vtkLinearExtrusionFilter

from vtkbool.vtkBool import vtkPolyDataBooleanFilter

a = 4.3333
b = 6.5

e = .15
f = .625

seqs = [
    [(2, 2, 2, 2), (3, 4, 2, 3)],
    [(3, 2, 3), (2, 2, 2, 2, 2, 2)],
    [(2, 2, 2, 2), (3, 2, 4, 3)],
    [(3, 2, 3), (2, 2, 2, 2, 2, 2)]
]

counts_a = [ sum(s_a) for s_a, _ in seqs ]
counts_b = [ sum(s_b) for _, s_b in seqs ]

assert counts_a.count(counts_a[0]) == len(counts_a)
assert counts_b.count(counts_b[0]) == len(counts_b)

count = 15

def transform(line, ang, dx, dy):
    tranformed = []

    for p in line:
        x = math.cos(ang)*p[0]-math.sin(ang)*p[1]
        y = math.sin(ang)*p[0]+math.cos(ang)*p[1]

        tranformed.append((x+dx, y+dy))

    return tranformed

def get_line(w, seq):
    x = (w+e)/sum(seq)

    xs = [0]

    for c in seq[:-1]:
        xs.append(xs[-1]+c*x)

    line = []

    for s, l in zip(xs, map(lambda c: c*x, seq)):
        line.extend([(s, 0), (s+l-e, 0), (s+l-e, e), (s+l, e)])

    del line[-3:]

    return line

def get_bricks_line(a, seq_a, b, seq_b):
    line_a = get_line(a, seq_a)
    line_b = get_line(b, seq_b)

    pl = line_a + transform(line_b, math.pi/2, a, 0) + transform(line_a, math.pi, a, b) + transform(line_b, 3*math.pi/2, 0, b)

    print(pl)

    return pl

def extrude_from_line(line, h, z=0):
    cell = vtkIdList()

    pts = vtkPoints()
    pts.SetDataTypeToDouble()

    [ (pts.InsertNextPoint(pt[0], pt[1], z), cell.InsertNextId(i)) for i, pt in enumerate(line) ]

    pd = vtkPolyData()
    pd.Allocate(1)
    pd.SetPoints(pts)
    pd.InsertNextCell(VTK_POLYGON, cell)

    extr = vtkLinearExtrusionFilter()
    extr.SetInputData(pd)
    extr.SetVector(0, 0, h)

    normals = vtkPolyDataNormals()
    normals.SetInputConnection(extr.GetOutputPort())
    normals.AutoOrientNormalsOn()

    return normals

def union(f_a, f_b):
    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, f_a.GetOutputPort())
    bf.SetInputConnection(1, f_b.GetOutputPort())

    return bf

all_seqs = []

while len(all_seqs) < count:
    all_seqs.extend(seqs)

    del all_seqs[count:]

bfs = []

for i, (s_a, s_b) in enumerate(all_seqs):
    bricks_line = get_bricks_line(a, s_a, b, s_b)

    bricks = extrude_from_line(bricks_line, f-e, i*f)
    gap = extrude_from_line([(e, e), (a-e, e), (a-e, b-e), (e, b-e)], e, (i+1)*f-e)

    bfs.append(union(bricks, gap))

bf = reduce(lambda a, b: union(a, b), bfs)

clean = vtkCleanPolyData()
clean.SetInputConnection(bf.GetOutputPort())

chimney = union(clean, extrude_from_line([(0, 0), (a, 0), (a, b), (0, b)], .5, count*f))

writer = vtkPolyDataWriter()
writer.SetInputConnection(chimney.GetOutputPort())
writer.SetFileName('bricks.vtk')
writer.Update()

stl_writer = vtkSTLWriter()
stl_writer.SetInputConnection(chimney.GetOutputPort())
stl_writer.SetFileName('bricks.stl')
stl_writer.Update()
