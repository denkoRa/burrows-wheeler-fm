def bwt(s, sa):
    """Given string s returns BWT(s), using suffix array"""
    bw = []
    for idx in sa:
        if idx == 0:
            bw.append("$")
        else:
            bw.append(s[idx - 1])
    return "".join(bw)


    
    