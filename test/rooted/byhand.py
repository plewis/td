from math import pow

#         ref  other
kf  = pow( 1  -  3   , 2) # A|BCDE
kf += pow( 1  -  3   , 2) # B|ACDE
kf += pow( 3  -  6   , 2) # C|ABDE
kf += pow( 5  -  4   , 2) # D|ABCE
kf += pow( 5  -  4   , 2) # E|ABCD
kf += pow( 2  -  6   , 2) # AB|CDE
kf += pow( 4  -  0   , 2) # ABC|DE
kf += pow( 2  -  2   , 2) # DE|ABC
kf += pow( 0  -  3   , 2) # CDE|AB
print('KF = %.6f' % kf)

rf  = pow( 1  -  0   , 2) # ABC|DE
rf += pow( 0  -  1   , 2) # CDE|AB
print('RF = %.6f' % rf)
