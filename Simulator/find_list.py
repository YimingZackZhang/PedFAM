import sys

with open('sim_Parent.dat', 'r') as f:
    data = []
    for x in f.readlines():
        x=x.strip('\n')
        x=x.split(' ')
        data.append(x)

indv = int(sys.argv[1])
posize = int(sys.argv[2])
gen = [indv-1]
for i in range(1,len(data)):
    
    gen_1 = []
    for j in range(len(gen)):
        gen_1.append(int(data[len(data)-i-1][gen[j]*2]))
        gen_1.append(int(data[len(data)-i-1][gen[j]*2 + 1]))
    gen = gen_1[:]
    #print(gen)
xx = gen[:]

a = 0
b = 0
xs = []
for i in range(len(xx)):
    if int(xx[i]+1) >= posize+1:
        b += 1
        xs.append('1')
    else:
        a += 1
        xs.append('0')
print('Number of founders from population A:',a,'Number of founders from population B:',b)
print('The ground truth configuration:',xs)



