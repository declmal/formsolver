import numpy as np


class Node:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
    def get_coord(self):
        return self.x, self.y, self.z
    def serialize(self):
        return ','.join([str(self.x), str((self.y)), str(self.z)])


def yaw(theta):
    cs = np.cos(theta)
    si = np.sin(theta)
    return np.array([
        [cs,-si,0],
        [si,cs,0],
        [0,0,1],
    ])

def pitch(theta):
    cs = np.cos(theta)
    si = np.sin(theta)
    return np.array([
        [cs,0,si],
        [0,1,0],
        [-si,0,cs],
    ])

def roll(theta):
    cs = np.cos(theta)
    si = np.sin(theta)
    return np.array([
        [1,0,0],
        [0,cs,-si],
        [0,si,cs],
    ])

def rotate_3d_vec(vec, yawMat, pitchMat, rollMat):
    return np.matmul(
        yawMat, np.matmul(
            pitchMat, np.matmul(rollMat, vec)
        )
    )


class Element:
    def __init__(self, nodes):
        self.nodes = nodes
    def get_nodes(self):
        return self.nodes
    def serialize(self, otype="default"):
        if otype == "default":
            return '\n'.join(
                [node.serialize() for node in self.nodes]
            )
        elif otype == "form":
            return '\n'.join(
                [node.serialize()+',' for node in self.nodes]
            )
        elif otype == "inp":
            return '\n'.join(
                [
                    str(i+1)+','+node.serialize()
                    for i, node in enumerate(self.nodes)
                ]
            )
        else:
            raise NotImplementedError
    @staticmethod
    def rotate_3d(elem, origin=Node(0,0,0), alpha=0.0, beta=0.0, gamma=0.0):
        nnodes = []
        xo, yo, zo = origin.get_coord()
        yawMat = yaw(alpha)
        pitchMat = pitch(beta)
        rollMat = roll(gamma)
        for node in elem.get_nodes():
            x, y, z = node.get_coord()
            vec = np.array([x-xo,y-yo,z-zo])
            nvec = rotate_3d_vec(vec, yawMat, pitchMat, rollMat)
            ndx, ndy, ndz = nvec.tolist()
            nx, ny, nz = xo+ndx, yo+ndy, zo+ndz
            nnodes.append(Node(nx,ny,nz))
        return type(elem)(nnodes)

# C3D20
c3d20rUnitNodes = [
    Node(x,y,z) for x, y, z in [
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1],
        [0.5, 0, 0],
        [1, 0.5, 0],
        [0.5, 1, 0],
        [0, 0.5, 0],
        [0.5, 0, 1],
        [1, 0.5, 1],
        [0.5, 1, 1],
        [0, 0.5, 1],
        [0, 0, 0.5],
        [1, 0, 0.5],
        [1, 1, 0.5],
        [0, 1, 0.5],
    ]
]
c3d20Unit = Element(c3d20rUnitNodes)

# C3D8
c3d8UnitNodes = [
    Node(x,y,z) for x,y,z in [
      [0, 0, 0],
      [1, 0, 0],
      [1, 1, 0],
      [0, 1, 0],
      [0, 0, 1],
      [1, 0, 1],
      [1, 1, 1],
      [0, 1, 1],
    ]
]
c3d8Unit = Element(c3d8UnitNodes)

if __name__ == '__main__':
    origin = Node(0.5,0.5,0.5)
    alpha = 0.3
    beta = 1.5
    gamma = -2.1
    # elem = Element.rotate_3d(
        # c3d8Unit, origin=origin,
        # alpha=alpha, beta=beta, gamma=gamma)
    # print(elem.serialize())
    otype = "inp"
    elem = Element.rotate_3d(
        c3d20Unit, origin=origin,
        alpha=alpha, beta=beta, gamma=gamma)
    print(elem.serialize(otype=otype))
