s = 48
p2s = 2 ** s
di1 = [3, 7, 15, 21, 27, 35, 59, 65, 77, 87, 93]
di2 = [5, 9, 17, 23, 33, 45, 63, 75, 83, 89, 95]
B = []
for i in range(11):
    B.append(p2s - di1[i])
    B.append(p2s - di2[i])
print(B)