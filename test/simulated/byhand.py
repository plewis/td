from math import pow

#         true - smc
ss  = pow(.023484-.009367,2) # A|BCCDE
ss += pow(.010832-.005032,2) # B|ACDE
ss += pow(.054362-.033813,2) # C|ABDE
ss += pow(.010832-.005032,2) # D|ABCE
ss += pow(.025723-.009101,2) # E|ABCD
ss += pow(.012652-.004069,2) # BD|ACE
ss += pow(.002239-.000000,2) # ABD|CE
ss += pow(.000000-.000266,2) # BDE|AC
ss += pow(.028638-.024446,2) # ABDE|C
print('%.6f' % ss)
