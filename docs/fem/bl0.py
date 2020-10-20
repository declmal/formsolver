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
    BL0 = []
    for k in range(6):
        ijs = k2ij[k]
        row = []
        for a in range(3*m_, 3*m_+3):
            m = a // 3
            p = a - 3*m
            m = "m"
            elem = {}
            for i, j in ijs:
                if j != p:
                    continue
                key = "_0 h_{},{}".format(m, i)
                if key not in elem:
                    elem[key] = 1
                else:
                    elem[key] += 1
            if not elem:
                elem = "0"
            else:
                elem = ["{} {}".format(cnt, key) if cnt > 2 else key for key, cnt in elem.items()]
                elem = '  +  '.join(elem)
            row.append(elem)
        BL0.append(row)
    print("BL0\n")
    for row in BL0:
        nrow = ''.join(["{:>15}".format(elem) for elem in row])
        print(nrow)
        print("")

