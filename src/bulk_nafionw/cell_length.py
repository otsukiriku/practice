
filestep = 400

for i in range(0,100):
    file_name = 'dump.pos.'+str(i*filestep)
    print(file_name)
    path = './'+file_name
    lst=[]
    with open(path, mode='r') as f:
        for j in range (1,9):
            z_cell = f.readline().split()
        length=float(z_cell[1])-float(z_cell[0])
        print(length)
        lst.append(length)
    
    print(max(lst))
