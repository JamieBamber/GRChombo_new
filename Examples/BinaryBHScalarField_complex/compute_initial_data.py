import math
# 

## Compute the bare mass from the ADM masses
"""
the ADM masses of the black holes are M1 M2
q = M2/M1 where M2 >= M1
M = M1 + M2
pt = tangential momentum
pr = radial momentum
D = distance between the black holes (in units of M)

m1 = M1 - P1^2/(8M1) - M1*M2/(2*D)
"""

M = 1
q = 2
D = 10
pt = 0.085599
pr = 0.0007948
p2 = pt**2 + pr**2

M1 = M/(1+q)
M2 = q*M/(1+q)

m1 = M1 - p2/(8*M1) - M1*M2/(2*D)
m2 = M2 - p2/(8*M2) - M1*M2/(2*D)

print("m1 = ", m1)
print("m2 = ", m2)
