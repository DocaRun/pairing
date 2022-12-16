"""
create : 2022/12/1
@author: DocaRun
"""

import David_reduction


# Prime
e1 = 0xfa
e2 = 0x9f
p = (2 ** e1) * (3 ** e2) - 1
e = 503

for i in range(10000):
    X, XR, outint = David_reduction.mainReduction()
    print('X',X % p)
    print('XR', XR)
    print('outint',outint)
    print('i', i)
    succ = 0
    pati = 0
    if (X % p == outint):
        succ += 1
    elif(X % p == outint - p):
        pati += 1
    else:
        exit()
    

print('correct reduction case', succ)
print('partial reduction case', pati)