import os
import math
import logging
import subprocess
import numpy as np
import utils

from os import path
from rtree import index
from utils import timethis


class Coord:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def get_coord(self):
        return self.x, self.y, self.z

    def get_coord_bbox(self):
        return self.get_coord() + self.get_coord()

    def __str__(self):
        return "({}, {}, {})".format(self.x, self.y, self.z)

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
        # print(coordBbox)
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
        assert len(nodeIds) == 20, len(nodeIds)
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
    def __init__(self, geoName, ElemType, **kwargs):
        attrs = kwargs.get('attrs')
        logger = logging.getLogger("form.grid.{}".format(geoName))
        for k, v in attrs.items():
            logger.info("attribute: {:>15}, value: {}".format(k, v))
        self.name = geoName
        self.ElemType = ElemType
        self.elemSet = ElemSet()
        self.nodeSet = NodeSet()


class TunnelLining(GeoType):
    @timethis
    def __init__(
        self, geoName, ElemType,
        depth=100000, num_layers=500,
        r_out=9000, r_in=8500, num_loops=3, num_slices=50,
        org_coord=Coord(0,0,0), **kwargs):
        super().__init__(geoName, ElemType, **kwargs)
        # create elems
        dy = depth / num_layers
        dr = (r_out-r_in) / num_loops
        dphi = np.pi*2 / num_slices
        ox, oy, oz = org_coord.get_coord()
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
        # check elems
        numNodesRef = 0
        ElemTypeName = ElemType.get_type()
        if ElemTypeName.startswith("C3D8"):
            numNodesRef += (num_layers+1) * (num_loops+1) * num_slices
        elif ElemTypeName.startswith("C3D20"):
            numNodesRef += (num_layers+1) * (num_loops+1) * num_slices * 2
            numNodesRef += (num_layers+1) * num_loops * num_slices
            numNodesRef += num_layers * (num_loops+1) * num_slices
        numNodes = self.nodeSet.get_size()
        assert numNodes == numNodesRef, \
            "num_nodes: {}, numNodesRef: {}".format(numNodes, numNodesRef)

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


class QuadVtk:
    def __init__(self, coordIds):
        self.coordIds = coordIds

    def get_coordIds(self):
        return self.coordIds


class SurroundingRock(GeoType):
    @timethis
    def __init__(
        self, geoName, ElemType,
        depth=100000, num_layers=500,
        radius=8500, num_elem_arc=45,
        side_length=17000, num_elem_side=24,
        l_top=8500, num_elem_top=45,
        l_right=8000, num_elem_right=45,
        l_bottom=7500, num_elem_bottom=40,
        l_left=7000, num_elem_left=40,
        org_coord=Coord(0,0,0), **kwargs):
        super().__init__(geoName, ElemType, **kwargs)
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
        # gmsh to vtk
        ext = 'vtk'
        outPath = path.join(gmshDir, geoName+'.'+ext)
        runstr = "gmsh {} -2 -o {} -format {}".format(geoPath, outPath, ext)
        subprocess.check_call(runstr, timeout=20, shell=True)
        # parse vtk file and check validity
        self.nodeSet = NodeSet()
        num_points = 0
        with open(outPath, 'r') as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if line.startswith("POINTS"):
                num_points = eval(line.split(' ')[1])
                break
        assert num_points > 0
        coords_2d = []
        for j in range(num_points):
            i += 1
            x, y, z = [eval(v) for v in lines[i].split(' ')]
            coord = Coord(x, y, z)
            coords_2d.append(coord)
        num_cells = 0
        for j in range(i, len(lines)):
            if lines[j].startswith("CELLS"):
                num_cells = eval(lines[j].split(' ')[1])
                break
        assert num_cells > 0
        quads = []
        for i in range(num_cells):
            j += 1
            if lines[j].startswith("4") or lines[j].startswith("8"):
                quad = QuadVtk([eval(v) for v in lines[j].split(' ')][1:5])
                quads.append(quad)
        for i in range(j, len(lines)):
            if lines[i].startswith("CELL_TYPES"):
                num_cell_types = eval(lines[i].split(' ')[1])
                break
        assert num_cells == num_cell_types
        quad_types = []
        num_mid_points = 0
        for j in range(num_cell_types):
            i += 1
            quad_type = eval(lines[i])
            if quad_type == 21:
                num_mid_points += 1
            if quad_type not in [1, 3, 21]:
                quad_types.append(quad_type)
        assert len(quad_types) == len(quads) and \
            all([quad_type == 9 for quad_type in quad_types]) or \
            all([quad_type == 23 for quad_type in quad_types])
        xs, zs, xs0, zs0, ids = set(), set(), set(), set(), set()
        for quad in quads:
            for cid in quad.get_coordIds():
                ids.add(cid)
                x, z, _ = coords_2d[cid].get_coord()
                if x == side_length:
                    zs.add(z)
                if z == side_length:
                    xs.add(x)
                if x == 0:
                    zs0.add(z)
                if z == 0:
                    xs0.add(x)
        assert len(xs) == len(zs) == num_elem_side+1, \
            "len(xs): {}, len(zs): {}, num_elem_side: {}".format(
                len(xs), len(zs), num_elem_side)
        assert len(xs0) == len(zs0) and len(xs0) > 1
        #  xs, zs = sorted(list(xs)), sorted(list(zs))
        # assert len(xs) > 1 and len(zs) > 1
        # check_xs = [self._check_equal(xs[i+1]-xs[i], ds) \
            # for i in range(len(xs)-1)]
        # check_zs = [self._check_equal(zs[i+1]-zs[i], ds) \
            # for i in range(len(zs)-1)]
        # assert all(check_xs) and all(check_zs), \
            # "len xs: {}, len zs: {}, xs: {}, zs: {}, ds: {}".format(
                # len(xs), len(zs), xs, zs, ds)
        ds = side_length / num_elem_side
        ds0 = (side_length-radius) / (len(xs0)-1)
        for quad in quads:
            for cid in quad.get_coordIds():
                x, z, _ = coords_2d[cid].get_coord()
                if x == side_length:
                    coords_2d[cid] = Coord(x, round(z/ds)*ds, 0)
                if z == side_length:
                    coords_2d[cid] = Coord(round(x/ds)*ds, z, 0)
                if x == 0:
                    coords_2d[cid] = Coord(x, radius+round((z-radius)/ds0)*ds0, 0)
                if z == 0:
                    coords_2d[cid] = Coord(radius+round((x-radius)/ds0)*ds0, z, 0)
        # create elems
        dy = depth / num_layers
        ox, oy, oz = org_coord.get_coord()
        cy = oy
        dphi = np.pi / 2
        ls = [l_top, l_right, l_bottom, l_left]
        num_elems = [
            num_elem_top, num_elem_right, num_elem_bottom, num_elem_left]
        for i in range(num_layers):
            ny = cy + dy
            cphi = 0
            for j in range(4):
                # div j-0
                for quad in quads:
                    coords = []
                    for y in [cy, ny]:
                        for cid in quad.get_coordIds():
                            dx, dz, _ = coords_2d[cid].get_coord()
                            coord = self._get_coord(dx, dz, ox, oz, cphi, y)
                            coords.append(coord)
                    self.elemSet.create_elem(ElemType, self.nodeSet, coords)
                # div j-1, j-2, j-3
                for lx, num_x, lz, num_z, x, z in [
                    (
                        side_length,
                        num_elem_side,
                        ls[j],
                        num_elems[j],
                        0,
                        side_length
                    ),
                    (
                        ls[(j+1)%4],
                        num_elems[(j+1)%4],
                        side_length,
                        num_elem_side,
                        side_length,
                        0
                    ),
                    (
                        ls[(j+1)%4],
                        num_elems[(j+1)%4],
                        ls[j],
                        num_elems[j],
                        side_length,
                        side_length,
                    ),
                ]:
                    dx = lx / num_x
                    dz = lz / num_z
                    for k in range(num_x):
                        for l in range(num_z):
                            coords = [
                                self._get_coord(
                                    x+k*dx, z+l*dz, ox, oz, cphi, cy),
                                self._get_coord(
                                    x+(k+1)*dx, z+l*dz, ox, oz, cphi, cy),
                                self._get_coord(
                                    x+(k+1)*dx, z+(l+1)*dz, ox, oz, cphi, cy),
                                self._get_coord(
                                    x+k*dx, z+(l+1)*dz, ox, oz, cphi, cy),
                                self._get_coord(
                                    x+k*dx, z+l*dz, ox, oz, cphi, ny),
                                self._get_coord(
                                    x+(k+1)*dx, z+l*dz, ox, oz, cphi, ny),
                                self._get_coord(
                                    x+(k+1)*dx, z+(l+1)*dz, ox, oz, cphi, ny),
                                self._get_coord(
                                    x+k*dx, z+(l+1)*dz, ox, oz, cphi, ny),
                            ]
                            self.elemSet.create_elem(
                                ElemType, self.nodeSet, coords)
                cphi += dphi
            cy = ny
        # check elems
        numNodesRef = 0
        numEndNodesPlane = len(ids) * 4
        for j in range(4):
            for num_x, num_z in [
                (num_elem_side, num_elems[j]),
                (num_elems[(j+1)%4], num_elem_side),
                (num_elems[(j+1)%4], num_elems[j]),
            ]:
                numEndNodesPlane += (num_x+1) * (num_z+1)
        for j in range(4):
            numEndNodesPlane -= \
                2*(num_elems[j]+1) + (num_elems[(j+1)%4]+1) + \
                2*(num_elem_side+1) + len(xs0)
        numEndNodesPlane += 8
        ElemTypeName = ElemType.get_type()
        if ElemTypeName.startswith("C3D8"):
            numNodesRef += numEndNodesPlane * (num_layers+1)
        elif ElemTypeName.startswith("C3D20"):
            numMidNodesPlane = num_mid_points * 4
            for j in range(4):
                for num_x, num_z in [
                    (num_elem_side, num_elems[j]),
                    (num_elems[(j+1)%4], num_elem_side),
                    (num_elems[(j+1)%4], num_elems[j]),
                ]:
                    numMidNodesPlane += (num_x+1)*num_z + num_x*(num_z+1)
            for j in range(4):
                numMidNodesPlane -= \
                    2*num_elems[j] + num_elems[(j+1)%4] + \
                    2*num_elem_side + len(xs0) - 1
            numNodesRef += numEndNodesPlane * (2*num_layers+1)
            numNodesRef += numMidNodesPlane * (num_layers+1)
        numNodes = self.nodeSet.get_size()
        assert numNodes == numNodesRef, \
            "num_nodes: {}, numNodesRef: {}".format(numNodes, numNodesRef)

    def _check_equal(self, v, ref_v, tol=4e-6):
        return abs(v-ref_v) <= tol

    def _get_coord(self, dx, dz, ox, oz, phi, y):
        cs = np.cos(phi)
        si = np.sin(phi)
        rot_mat = np.array([[cs,si],[-si,cs]])
        v = np.array([dx, dz])
        nv = np.matmul(rot_mat, v)
        ndx, ndz = nv.tolist()
        return Coord(ox+ndx, y, oz+ndz)


class Model:
    def __init__(self):
        self.geoLst = []

    def add_geo(self, geo):
        self.geoLst.append(geo)

    @timethis
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

    @timethis
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
    # log init
    utils.log_init()
    logger = logging.getLogger("form.grid")
    # geometric attributes
    model = Model()
    depth = 1200
    r_in = 8500
    r_out = 9000
    esize_bd = 500
    esize_r_out = 350
    ElemType = C3D20
    org_coord = Coord(0,0,0)
    # surrounding rock
    geoName = "SurroundingRock"
    side_length = 1.3 * r_out
    l_top = 8500
    l_right = 8000
    l_bottom = 7500
    l_left = 6000
    esize_layer_rock = 400
    geoAttrs = {
        'depth': depth,
        'num_layers': round(depth/esize_layer_rock),
        'radius': r_out,
        'num_elem_arc': round(np.pi/2*r_out/esize_r_out),
        'side_length': side_length,
        'num_elem_side': round(0.5*side_length/esize_bd)*2,
        'l_top': l_top,
        'num_elem_top': round(l_top/esize_bd),
        'l_right': l_right,
        'num_elem_right': round(l_right/esize_bd),
        'l_bottom': l_bottom,
        'num_elem_bottom': round(l_bottom/esize_bd),
        'l_left': l_left,
        'num_elem_left': round(l_left/esize_bd),
        'org_coord': org_coord,
    }
    geo = SurroundingRock(geoName, ElemType, **geoAttrs, attrs=geoAttrs)
    model.add_geo(geo)
    # tunnel lining
    geoName = "TunnelLining"
    esize_slice = 160
    esize_loop = 160
    esize_layer_lining = 160
    geoAttrs = {
        'depth': depth,
        'num_layers': round(depth/esize_layer_lining),
        'r_out': r_out,
        'r_in': r_in,
        'num_loops': round((r_out-r_in)/esize_slice),
        'num_slices': round(0.25*np.pi*(r_in+r_out)/esize_loop) * 4,
        'org_coord': org_coord,
    }
    geo = TunnelLining(geoName, ElemType, **geoAttrs, attrs=geoAttrs)
    model.add_geo(geo)
    # output
    model.to_vtk()
