if __name__ == '__main__':
    S = []
    for k in range(9):
        SS = []
        # k = 3r + i
        r = k // 3
        i = k - 3*r
        for l in range(9):
            # l = 3s + j
            s = l // 3
            j = l - 3*s
            if r != s:
                SS.append("0")
            else:
                SS.append("S{}{}".format(i, j))
        S.append(SS)
    print("S0\n")
    for SS in S:
        nrow = ''.join(["{:>5}".format(elem) for elem in SS])
        print(nrow)
        print("")
