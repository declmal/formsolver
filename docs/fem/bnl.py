if __name__ == '__main__':
    m_ = 8
    start_ = 0
    S = []
    for k in range(9):
        SS = []
        r = k // 3
        i = k - 3*r
        for a in range(3*m_):
            m = a // 3
            p = a - 3*m
            if p != r:
                SS.append("0")
            else:
                SS.append("h{}{}".format(m, i))
        S.append(SS)
    print("BNL\n")
    for SS in S:
        nrow = ''.join(["{:>7}".format(elem) for elem in SS])
        print(nrow)
        print("")
