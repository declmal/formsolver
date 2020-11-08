import os
import math
import numpy as np

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

    @staticmethod
    def get_type():
        return "C3D20"

    @staticmethod
    def to_elem(nodeSet, coords, elemId):
        ncoords = [
            Coord.get_mid_coord(coords[0],coords[1]),
            Coord.get_mid_coord(coords[1],coords[2]),
            Coord.get_mid_coord(coords[2],coords[3]),
            Coord.get_mid_coord(coords[3],coords[1]),
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

    def create_elem(self, elemType, nodeSet, coords):
        elemId = len(self.elemLst)
        elem = elemType.to_elem(nodeSet, coords, elemId)
        self.elemLst.append(elem)

    def __iter__(self):
        return iter(self.elemLst)


class GeoType:
    def create_mesh(self):
        raise NotImplementedError


class TunnelGeoType(GeoType):
    def __init__(
        self, depth=100000, num_layers=500, r_out=9000, r_in=8500,
        num_loops=3, num_slices=50, orgCoord=Coord(0,0,0)):
        print("num element:", num_layers*num_loops*num_slices)
        self.depth = depth
        self.num_layers = num_layers
        self.r_out = r_out
        self.r_in = r_in
        self.num_loops = num_loops
        self.num_slices = num_slices
        self.orgCoord = orgCoord

    def create_mesh(self, elemType):
        dy = self.depth / self.num_layers
        dr = (self.r_out-self.r_in) / self.num_loops
        dphi = 360 / self.num_slices
        ox, oy, oz = self.orgCoord.get_coord()
        cy = oy
        elemSet = ElemSet()
        nodeSet = NodeSet()
        for i in range(self.num_layers):
            ny = cy + dy
            cr = self.r_in
            for j in range(self.num_loops):
                nr = cr + dr
                cphi = 0
                for k in range(self.num_slices):
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
                    elemSet.create_elem(elemType, nodeSet, coords)
                    cphi = nphi
                cr = nr
            cy = ny
        return elemSet, nodeSet

    def _get_coord(self, y, r, phi):
        x = r * np.sin(phi)
        z = r * np.cos(phi)
        return Coord(x, y, z)


class Model:
    def __init__(self, geo, elemType):
        self.elemType = elemType
        self.elemSet, self.nodeSet = geo.create_mesh(elemType)

    def to_inp(self, modelName="model", modelDir=path.expanduser("~/.model")):
        inpContent = "*Node\n"
        inpContent += '\n'.join(
            [node.to_inp() for node in self.nodeSet]
        )
        inpContent += "\n" + ','.join([
            "*Element","TYPE={}".format(elemType.get_type())
        ]) + "\n"
        inpContent += '\n'.join(
            [elem.to_inp() for elem in self.elemSet]
        )
        modelPath = path.join(modelDir, modelName+'.inp')
        with open(modelPath, 'w') as f:
            f.write(inpContent)

    def to_json(self, modelName="model", modelDir=path.expanduser("~/.model")):
        pass

if __name__ == "__main__":
    tunnelGeoAttrs = {
        'num_layers': 50,
        'num_loops': 3,
        'num_slices': 5,
    }
    geo = TunnelGeoType(**tunnelGeoAttrs)
    elemType = C3D20
    model = Model(geo, elemType)
    model.to_inp()
