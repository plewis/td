from math import pow

#         ref  other
ss  = pow( 1  -  3   , 2) # A|BCDE
ss += pow( 1  -  3   , 2) # B|ACDE
ss += pow( 3  -  6   , 2) # C|ABDE
ss += pow( 5  -  4   , 2) # D|ABCE
ss += pow( 5  -  4   , 2) # E|ABCD
ss += pow( 2  -  6   , 2) # AB|CDE
ss += pow( 4  -  0   , 2) # ABC|DE
ss += pow( 2  -  2   , 2) # DE|ABC
ss += pow( 0  -  3   , 2) # CDE|AB
print('%.6f' % ss)
