if __name__ == '__main__':
    kf = {
        0: 1,
        1: 1,
        2: 1,
        3: 2,
        4: 2,
        5: 2,
    }
    k2ij = {
        0: [(0,0)],
        1: [(1,1)],
        2: [(2,2)],
        3: [(0,1)],
        4: [(1,2)],
        5: [(0,2)],
    }
    m_ = 3
    BL0 = []
    for k in range(6):
        ijs = k2ij[k]
        factor = kf[k]
        row = []
        for a in range(3*m_, 3*m_+3):
            m = a // 3
            p = a - 3*m
            m = "m"
            elem = {}
            for i, j in ijs:
                keys = []
                if i == p:
                    key1 = "_0 h_{},{}".format(m, j)
                    keys.append(key1)
                if j == p:
                    key2 = "_0 h_{},{}".format(m, i)
                    keys.append(key2)
                for key in keys:
                    if key not in elem:
                        elem[key] = 1
                    else:
                        elem[key] += 1
            elem = {key: cnt*factor*0.5 for key, cnt in elem.items()}
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

