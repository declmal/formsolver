if __name__ == "__main__":
    dct = {
        0: [(0,0)],
        1: [(1,1)],
        2: [(2,2)],
        3: [(0,1),(1,0)],
        4: [(1,2),(2,1)],
        5: [(0,2),(2,0)],
    }
    mat = []
    for k in range(6):
        row = []
        for l in range(6):
            cnt_lambda = 0
            cnt_mu = 0
            for i, j in dct[k]:
                for r, s in dct[l]:
                    if i == j and r == s:
                        cnt_lambda += 1
                    if i == r and j == s:
                        cnt_mu += 1
                    if i == s and j == r:
                        cnt_mu += 1
            items = []
            if cnt_lambda == 1:
                items.append("L")
            elif cnt_lambda > 1:
                items.append(str(cnt_lambda) + "L")
            if cnt_mu == 1:
                items.append("M")
            elif cnt_mu > 1:
                items.append(str(cnt_mu) + "M")
            if items:
                row.append('+'.join(items))
            else:
                row.append('0')
        mat.append(row)
    for row in mat:
        print(','.join(["{:>10}".format(e) for e in row]))
