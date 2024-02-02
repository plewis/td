from math import pow

#          ref  other
ss  = pow(  1 -   3   , 2) # A|BCDE
ss += pow(  1 -   3   , 2) # B|ACDE
ss += pow(  3 -   6   , 2) # C|ABDE
ss += pow(  5 -   4   , 2) # D|ABCE
ss += pow(  5 -   4   , 2) # E|ABCD
ss += pow(  2 -   9   , 2) # AB|CDE
ss += pow(  6 -   2   , 2) # ABC|DE
print('%.6f' % ss)
