if __name__ == '__main__':
    k2ij = {
        0: [(0,0)],
        1: [(1,1)],
        2: [(2,2)],
        3: [(0,1), (1,0)],
        4: [(1,2), (2,1)],
        5: [(0,2), (2,0)],
    }
    m_ = 3
    BL1 = []
    for k in range(6):
        ijs = k2ij[k]
        row = []
        for a in range(3*m_, 3*m_+3):
            m = a // 3
            p = a - 3*m
            m = "m"
            # m = 0
            elem = {}
            for i, j in ijs:
                key = "_0 h_{},{}  _0^t u_{},{}".format(m, i, p, j)
                # key = "_0^t u_{},{} _0 h_{},{}".format(p+1, j+1, m+1, i+1)
                if key not in elem:
                    elem[key] = 1
                else:
                    elem[key] += 1
            if not elem:
                elem = "0"
            else:
                elem = ["{} {}".format(cnt, key) if cnt > 2 else key for key, cnt in elem.items()]
                elem = '  +  '.join(elem)
                # elem = '  +  '.join(elem[::-1])
            row.append(elem)
        BL1.append(row)
    print("BL1\n")
    for row in BL1:
        nrow = ''.join(["{:>60}".format(elem) for elem in row])
        print(nrow)
        print("")

