from checkpoints import Checkpoints
from bw import bwt
from sa import suffix_array_basic, suffix_array_fast

class FMIndex:
    def cut_suffix_array(self, full_sa, sa_step):
        """Sizes down suffix array taking every sa_step-th element, the rest are removed"""
        res_sa = {}
        for i in range(len(full_sa)):
            if full_sa[i] % sa_step == 0:
                res_sa[i] = full_sa[i]
        return res_sa
    
    def calc_first_col(self, bwt):
        """
        Returns first column of BWT matrix.
        It is enough to keep only number of characters less than current, because it is sorted column.
        """
        cnts = {}
        for c in bwt:
            cnts[c] = cnts.get(c, 0) + 1
        ret = {}
        total = 0
        for c, cnt in sorted(cnts.items()):
            ret[c] = total
            total += cnt
        return ret
    
    def __init__(self, T, suffix_array, cp_step, sa_step):
        if T[-1] != '$':
            T += '$'
        full_sa = suffix_array(T)
        self.bwt = bwt(T, full_sa)
        self.sa = self.cut_suffix_array(full_sa, sa_step)
        self.cps = Checkpoints(self.bwt, cp_step)
        self.first_col = self.calc_first_col(self.bwt)
        
    def count(self, c):
        """Return number of characters less than c"""
        if c not in self.first_col:
            for cc in sorted(self.first_col.keys()):
                if c < cc: 
                    return self.first_col[cc]
            return self.first_col[cc]
        else:
            return self.first_col[c]
        
    def interval(self, p):
        """Returns inclusive interval of BWM rows where p is prefix"""
        l, r = 0, len(self.bwt) - 1
        print (self.bwt)
        print (self.first_col)
        for i in range(len(p) - 1, -1, -1): #start from last character and go backwards
            l = self.cps.rank(self.bwt, p[i], l - 1) + self.count(p[i])
            r = self.cps.rank(self.bwt, p[i], r) + self.count(p[i]) - 1
            if r < l:
                break
        return l, r + 1
    
    def get_offset(self, r):
        """Given row in matrix returns its real offset in text"""
        steps = 0
        while r not in self.sa:
            c = self.bwt[r]
            r = self.cps.rank(self.bwt, c, r - 1) + self.count(c)
            steps += 1
        return self.sa[r] + steps
    
    def find_occurences(self, p):
        """Returns all occurences of pattern p"""
        l, r = self.interval(p)
        return [self.get_offset(i) for i in range(l, r)]


#testing
test = "abaaba"
fm_index = FMIndex(test, suffix_array_fast, 1, 1)
fm_index.interval('ba')