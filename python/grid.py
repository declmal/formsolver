import os
import math
import numpy as np
import subprocess

from os import path
from rtree import index


class Coord:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def get_coord(self):
        return self.x, self.y, self.z

    def get_coord_bbox(self):
        return self.get_coord() + self.get_coord()

    @staticmethod
    def get_mid_coord(coord1, coord2):
        x1, y1, z1 = coord1.get_coord()
        x2, y2, z2 = coord2.get_coord()
        return Coord(0.5*(x1+x2), 0.5*(y1+y2), 0.5*(z1+z2))


class Node:
    def __init__(self, coord, nodeId):
        self.coord = coord
        self.nodeId = nodeId

    def get_coord(self):
        return self.coord.get_coord()

    def get_id(self):
        return self.nodeId

    def to_inp(self):
        x, y, z = self.get_coord()
        return ','.join([str(self.nodeId+1), str(x), str(y), str(z)])

    def to_vtk_point(self):
        return ' '.join([str(c) for c in self.get_coord()])

    def to_vtk_cell(self):
        return ' '.join([str(1), str(self.nodeId)])

    def to_vtk_cell_type(self):
        return str(1)


class NodeSet:
    def __init__(self, dstTol=1e-6):
        self.nodeLst = []
        self.dstTol = dstTol
        p = index.Property()
        p.dimension = 3
        self.idx = index.Index(properties=p)

    def insert_node(self, coord):
        cx, cy, cz = coord.get_coord()
        coordBbox = coord.get_coord_bbox()
        hits = list(self.idx.nearest(coordBbox, 1, objects=True))
        if hits:
            nx, ny, nz, _, _, _ = hits[0].bbox
            dx, dy, dz = nx-cx, ny-cy, nz-cz
            dst = math.sqrt(dx*dx + dy*dy + dz*dz)
            if dst < self.dstTol:
                return hits[0].id
        nodeId = len(self.nodeLst)
        node = Node(coord, nodeId)
        self.nodeLst.append(node)
        self.idx.insert(nodeId, coordBbox)
        print(coordBbox)
        return nodeId

    def get_size(self):
        return len(self.nodeLst)

    def __iter__(self):
        return iter(self.nodeLst)


class Element:
    def __init__(self, elemId, nodeIds):
        self.elemId = elemId
        self.nodeIds = nodeIds

    def to_inp(self):
        return ','.join(
            [str(self.elemId+1)] +
            [str(nodeId+1) for nodeId in self.nodeIds]
        )

    def to_vtk_cell(self):
        return ' '.join(
            [str(len(self.nodeIds))] + \
            [str(nodeId) for nodeId in self.nodeIds]
        )

    @classmethod
    def to_vtk_cell_type(cls):
        raise NotImplementedError

    @classmethod
    def get_num_points(cls):
        raise NotImplementedError

    @staticmethod
    def get_type():
        raise NotImplementedError

    @staticmethod
    def to_elem(nodeSet, coords, elemId):
        raise NotImplementedError


class C3D8(Element):
    def __init__(self, elemId, nodeIds):
        assert len(nodeIds) == 8
        super().__init__(elemId, nodeIds)

    @classmethod
    def to_vtk_cell_type(cls):
        return str(12)

    @classmethod
    def get_num_points(cls):
        return 8

    @staticmethod
    def get_type():
        return "C3D8"

    @staticmethod
    def to_elem(nodeSet, coords, elemId):
        return C3D8(elemId, [
            nodeSet.insert_node(coord) for coord in coords
        ])

class C3D8I(C3D8):
    @staticmethod
    def get_type():
        return "C3D8I"


class C3D20(Element):
    def __init__(self, elemId, nodeIds):
        assert len(nodeIds) == 20
        super().__init__(elemId, nodeIds)

    @classmethod
    def to_vtk_cell_type(cls):
        return str(25)

    @classmethod
    def get_num_points(cls):
        return 20

    @staticmethod
    def get_type():
        return "C3D20"

    @staticmethod
    def to_elem(nodeSet, coords, elemId):
        ncoords = [
            Coord.get_mid_coord(coords[0],coords[1]),
            Coord.get_mid_coord(coords[1],coords[2]),
            Coord.get_mid_coord(coords[2],coords[3]),
            Coord.get_mid_coord(coords[3],coords[0]),
            Coord.get_mid_coord(coords[4],coords[5]),
            Coord.get_mid_coord(coords[5],coords[6]),
            Coord.get_mid_coord(coords[6],coords[7]),
            Coord.get_mid_coord(coords[7],coords[4]),
            Coord.get_mid_coord(coords[0],coords[4]),
            Coord.get_mid_coord(coords[1],coords[5]),
            Coord.get_mid_coord(coords[2],coords[6]),
            Coord.get_mid_coord(coords[3],coords[7]),
        ]
        return C3D20(elemId, [
            nodeSet.insert_node(coord) for coord in coords+ncoords
        ])


class C3D20R(C3D20):
    @staticmethod
    def get_type():
        return "C3D20R"


class ElemSet:
    def __init__(self):
        self.elemLst = []

    def create_elem(self, ElemType, nodeSet, coords):
        elemId = len(self.elemLst)
        elem = ElemType.to_elem(nodeSet, coords, elemId)
        self.elemLst.append(elem)

    def get_size(self):
        return len(self.elemLst)

    def __iter__(self):
        return iter(self.elemLst)


class GeoType:
    def __init__(self, geoName, ElemType):
        self.name = geoName
        self.ElemType = ElemType
        self.elemSet = ElemSet()
        self.nodeSet = NodeSet()


class TunnelLining(GeoType):
    def __init__(
        self, geoName, ElemType,
        depth=100000, num_layers=500, r_out=9000, r_in=8500,
        num_loops=3, num_slices=50, orgCoord=Coord(0,0,0)):
        super().__init__(geoName, ElemType)
        print("geo name: {}, num element: {}".format(
            geoName, num_layers*num_loops*num_slices))
        dy = depth / num_layers
        dr = (r_out-r_in) / num_loops
        dphi = np.pi*2 / num_slices
        ox, oy, oz = orgCoord.get_coord()
        cy = oy
        for i in range(num_layers):
            ny = cy + dy
            cr = r_in
            for j in range(num_loops):
                nr = cr + dr
                cphi = 0
                for k in range(num_slices):
                    nphi = cphi + dphi
                    coords = [
                        self._get_coord(cy, cr, cphi),
                        self._get_coord(cy, cr, nphi),
                        self._get_coord(ny, cr, nphi),
                        self._get_coord(ny, cr, cphi),
                        self._get_coord(cy, nr, cphi),
                        self._get_coord(cy, nr, nphi),
                        self._get_coord(ny, nr, nphi),
                        self._get_coord(ny, nr, cphi),
                    ]
                    self.elemSet.create_elem(ElemType, self.nodeSet, coords)
                    cphi = nphi
                cr = nr
            cy = ny
        numNodes = self.nodeSet.get_size()
        ElemTypeName = ElemType.get_type()
        numNodesRef = 0
        if ElemTypeName.startswith("C3D8"):
            numNodesRef += (num_layers+1) * (num_loops+1) * num_slices
        elif ElemTypeName.startswith("C3D20"):
            numNodesRef += (num_layers+1) * (num_loops+1) * num_slices * 2
            numNodesRef += (num_layers+1) * num_loops * num_slices
            numNodesRef += num_layers * (num_loops+1) * num_slices
        assert numNodes == numNodesRef

    def _get_coord(self, y, r, phi):
        x = r * np.sin(phi)
        z = r * np.cos(phi)
        return Coord(x, y, z)

class GMSHItem:
    def __init__(self, gid):
        self.gid = gid

    def get_gid(self):
        return self.gid


class PointGMSH2D(Coord, GMSHItem):
    def __init__(self, gid, x, y, esize):
        Coord.__init__(self, x, y, 0)
        GMSHItem.__init__(self, gid)
        self.esize = esize

    def get_esize(self):
        return self.esize


class LineGMSH(GMSHItem):
    def __init__(self, gid, ipt1, ipt2):
        self.gid = gid
        self.ipt1 = ipt1
        self.ipt2 = ipt2

    def get_ipt(self):
        return self.ipt1, self.ipt2


class SurroundingRock(GeoType):
    def __init__(
        self, geoName, ElemType,
        radius=1, num_elem_arc=40, side_length=2, num_elem_side=24):
        super().__init__(geoName, ElemType)
        contents = []
        # points
        esize_arc = 2 * 0.5 * np.pi * radius / num_elem_arc
        esize_side = 2 * side_length / num_elem_side
        points = [
            PointGMSH2D(0, radius, 0, esize_arc),
            PointGMSH2D(1, 0, 0, esize_side),
            PointGMSH2D(2, 0, radius, esize_arc),
            PointGMSH2D(3, 0, side_length, esize_side),
            PointGMSH2D(4, side_length, side_length, esize_side),
            PointGMSH2D(5, side_length, 0, esize_side),
        ]
        for point in points:
            gid = point.get_gid()
            x, y, _ = point.get_coord()
            esize = point.get_esize()
            content = "Point({}) = {{{}, {}, 0, {}}};".format(
                gid, x, y, esize)
            contents.append(content)
        # circle
        contents.append("Circle(6) = {0, 1, 2};")
        # lines
        lines = [
            LineGMSH(7, 2, 3),
            LineGMSH(8, 3, 4),
            LineGMSH(9, 4, 5),
            LineGMSH(10, 5, 0),
        ]
        for line in lines:
            gid = line.get_gid()
            ipt1, ipt2 = line.get_ipt()
            content = "Line({}) = {{{}, {}}};".format(gid, ipt1, ipt2)
            contents.append(content)
        # line loop
        contents.append("Line Loop(11) = {6, 7, 8, 9, 10};")
        # plane surface
        contents.append("Plane Surface(12) = {11};")
        # gmsh configs
        contents.append('Mesh.RecombinationAlgorithm = 1; // blossom')
        contents.append('Mesh.RecombineAll = 1; // turns on quads')
        contents.append('Mesh.SubdivisionAlgorithm = 1; // quadrangles only')
        contents.append('Mesh.CharacteristicLengthExtendFromBoundary = 1;')
        contents.append('Mesh.CharacteristicLengthMin = 0;')
        contents.append('Mesh.CharacteristicLengthMax = 1e+022;')
        contents.append('Mesh.CharacteristicLengthFromPoints = 1;')
        contents.append('Mesh.Algorithm = 8; // delquad = delauny for quads')
        ElemTypeName = ElemType.get_type()
        if ElemTypeName.startswith("C3D8"):
            contents.append('Mesh.ElementOrder = 1; // linear or second set here')
        elif ElemTypeName.startswith("C3D20"):
            contents.append('Mesh.ElementOrder = 2; // linear or second set here')
            contents.append(
                'Mesh.SecondOrderIncomplete = 1; ' + \
                '// no face node w/ 2nd order')
        else:
            assert False
        contents.append('Mesh.SaveGroupsOfNodes = 1; // save node groups')
        # write to geo
        gmshDir = path.join(os.getcwd(), "build", "gmsh")
        os.makedirs(gmshDir, exist_ok=True)
        geoPath = path.join(gmshDir, geoName+".geo")
        with open(geoPath, 'w') as f:
            f.write('\n'.join(contents))
        # gmsh
        ext = 'vtk'
        outPath = path.join(gmshDir, geoName+'.'+ext)
        runstr = "gmsh {} -2 -o {} -format {}".format(geoPath, outPath, ext)
        subprocess.check_call(runstr, timeout=20, shell=True)
        # parse mesh file
        self.nodeSet = NodeSet()
        num_points = 0
        with open(outPath, 'r') as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if line.startswith("POINTS"):
                num_points = eval(line.split(' ')[1])
                break
        assert num_points > 0
        coords = []
        for j in range(num_points):
            i += 1
            x, y, z = [eval(v) for v in lines[i].split(' ')]
            coord = Coord(x, y, z)
            coords.append(coord)
        print(len(coords))
        exit()


class Model:
    def __init__(self):
        self.geoLst = []

    def add_geo(self, geo):
        self.geoLst.append(geo)

    def to_inp(self, modelName="model", modelDir=path.expanduser("~/.model")):
        for geo in self.geoLst:
            inpContent = "*Node\n"
            inpContent += '\n'.join(
                [node.to_inp() for node in geo.nodeSet]
            )
            inpContent += "\n" + ','.join([
                "*Element","TYPE={}".format(geo.ElemType.get_type())
            ]) + "\n"
            inpContent += '\n'.join(
                [elem.to_inp() for elem in geo.elemSet]
            )
            geoPath = path.join(modelDir, modelName+'-'+geo.name+'.inp')
            with open(geoPath, 'w') as f:
                f.write(inpContent)

    def to_vtk(self, modelName="model", modelDir=path.expanduser("~/.model")):
        for geo in self.geoLst:
            numElems = geo.elemSet.get_size()
            numNodes = geo.nodeSet.get_size()
            content = [
                "# vtk DataFile Version 2.0",
                modelName,
                "ASCII",
                "DATASET UNSTRUCTURED_GRID",
                ' '.join([
                    "POINTS", str(geo.nodeSet.get_size()), "double"
                ]),
                '\n'.join([
                    node.to_vtk_point() for node in geo.nodeSet
                ]) + '\n',
                ' '.join([
                    "CELLS",
                    str(numNodes+numElems),
                    str(
                        numNodes*2 + \
                        numElems*(1+geo.ElemType.get_num_points())
                    )
                ]),
                '\n'.join([
                    node.to_vtk_cell() for node in geo.nodeSet
                ]),
                '\n'.join([
                    elem.to_vtk_cell() for elem in geo.elemSet
                ]) + '\n',
                ' '.join([
                    "CELL_TYPES",
                    str(numNodes+numElems),
                ]),
                '\n'.join([
                    node.to_vtk_cell_type() for node in geo.nodeSet
                ]),
                '\n'.join([
                    elem.to_vtk_cell_type() for elem in geo.elemSet
                ]),
            ]
            geoPath = path.join(modelDir, modelName+'-'+geo.name+'.vtk')
            with open(geoPath, "w") as f:
                f.write('\n'.join(content))

    def to_json(self, modelName="model", modelDir=path.expanduser("~/.model")):
        pass

if __name__ == "__main__":
    model = Model()
    # surrounding rock
    geoName = "SurroundingRock"
    ElemType = C3D8
    geoAttrs = {
    }
    geo = SurroundingRock(geoName, ElemType, **geoAttrs)
    model.add_geo(geo)
    # tunnel lining
    geoName = "TunnelLining"
    ElemType = C3D8
    geoAttrs = {
        'depth': 1000,
        'num_layers': 5,
        'num_loops': 3,
        'num_slices': 350,
    }
    geo = TunnelLining(geoName, ElemType, **geoAttrs)
    model.add_geo(geo)
    # output
    model.to_vtk()
