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
    BL1 = []
    for k in range(6):
        ijs = k2ij[k]
        factor = kf[k]
        row = []
        for a in range(3*m_, 3*m_+3):
            m = a // 3
            p = a - 3*m
            m = "m"
            # m = 0
            elem = {}
            for i, j in ijs:
                key1 = "_0 h_{},{}  _0^t u_{},{}".format(m, j, p, i)
                key2 = "_0 h_{},{}  _0^t u_{},{}".format(m, i, p, j)
                for key in [key1, key2]:
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
                # elem = '  +  '.join(elem[::-1])
            row.append(elem)
        BL1.append(row)
    print("BL1\n")
    for row in BL1:
        nrow = ''.join(["{:>60}".format(elem) for elem in row])
        print(nrow)
        print("")

