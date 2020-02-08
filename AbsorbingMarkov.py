
#get the number of transient nodes based on the values!
def number_of_transients(matrix):
    if len(matrix) == 0:
        raise Exception("empty matrix")
    trs=[] #list of transient nodes(may conatin duplicate values)
    for r in range(len(matrix)):
        for c in range(len(matrix[r])):
            if matrix[r][c] != 0 and matrix[r][c] != 1:  #trasnient node
                trs.append(r)
    return len(trs)

#get the list of distinct transient nodes 
def getTransientNodes(matrix):
    if len(matrix) == 0:
        raise Exception("empty matrix")
    trs=[]
    for r in range(len(matrix)):
        for c in range(len(matrix[r])):
            if matrix[r][c] != 0 and matrix[r][c] != 1:
                trs.append(r)
    res=[]
    [res.append(x) for x in trs if x not in res] 
    return res

#get the number of absorbing nodes based on the values!
def number_of_absorbing(matrix):
    if len(matrix) == 0:
        raise Exception("empty matrix")
    abs=[]
    for r in range(len(matrix)):
        for c in range(len(matrix[r])):
            if matrix[r][c] == 1:
                abs.append(r)
    return len(abs)

#get the list of distinct absorbing nodes 
def getAbosrbingNodes(matrix):
    if len(matrix) == 0:
        raise Exception("empty matrix")
    abs=[]
    for r in range(len(matrix)):
        for c in range(len(matrix[r])):
            if matrix[r][c] == 1:
                abs.append(r)
    return abs

#produce the identity matrix with given size t
def identityMatrix(t):
    m = []
    for i in range(t):
        r = []
        for j in range(t):
            r.append(int(i == j))
        m.append(r)
    return m

#produce the zero matrix with given size t
def allzeroMatrix(t):
    m = []
    for i in range(t):
        r = []
        for j in range(t):
            r.append(0)
        m.append(r)
    return m

# check if the matrix is zero
def isZero(m):
    for r in range(len(m)):
        for c in range(len(m[r])):
            if m[r][c] != 0:
                return False
    return True

# subtract two matrices
def subtract(i, q):
    if len(i) != len(i[0]) or len(q) != len(q[0]):
        raise Exception("non-square matrices")

    if len(i) != len(q) or len(i[0]) != len(q[0]):
        raise Exception("Cannot subtract matrices of different sizes")
    s = []
    for r in range(len(i)):
        sRow = []
        for c in range(len(i[r])):
            sRow.append(i[r][c] - q[r][c])
        s.append(sRow)
    return s

# transpose matrix
def transposeMatrix(m):
    t = []
    for r in range(len(m)):
        tRow = []
        for c in range(len(m[r])):
            if c == r:
                tRow.append(m[r][c])
            else:
                tRow.append(m[c][r])
        t.append(tRow)
    return t

# multiply two matrices
def multiply(a, b):
    if a == [] or b == []:
        raise Exception("Cannot multiply empty matrices")

    if len(a[0]) != len(b):
        raise Exception("Cannot multiply matrices of incompatible sizes")

    m = []
    rows = len(a)
    cols = len(b[0])
    iters = len(a[0])

    for r in range(rows):
        mRow = []
        for c in range(cols):
            sum = 0
            for i in range(iters):
                sum += a[r][i]*b[i][c]
            mRow.append(sum)
        m.append(mRow)
    return m

def getMatrixMinor(m,i,j):
    return [row[:j] + row[j+1:] for row in (m[:i]+m[i+1:])]

# matrix determinant
def getMatrixDeternminant(m):
    #base case for 2x2 matrix
    if len(m) == 2:
        return m[0][0]*m[1][1]-m[0][1]*m[1][0]

    d = 0
    for c in range(len(m)):
        d += ((-1)**c)*m[0][c]*getMatrixDeternminant(getMatrixMinor(m,0,c))

    return d

#matrix Inverse
def getMatrixInverse(m):
    d = getMatrixDeternminant(m)

    if d == 0:
        raise Exception("Cannot get inverse of matrix with zero determinant")

    #special case for 2x2 matrix:
    if len(m) == 2:
        return [[m[1][1]/d, -1*m[0][1]/d],
                [-1*m[1][0]/d, m[0][0]/d]]

    #find matrix of cofactors
    cofactors = []
    for r in range(len(m)):
        cofactorRow = []
        for c in range(len(m)):
            minor = getMatrixMinor(m,r,c)
            cofactorRow.append(((-1)**(r+c)) * getMatrixDeternminant(minor))
        cofactors.append(cofactorRow)
    cofactors = transposeMatrix(cofactors)
    for r in range(len(cofactors)):
        for c in range(len(cofactors)):
            cofactors[r][c] = cofactors[r][c]/d
    return cofactors

#Rmatrix producer
def getRMatrix(trs,abso):
    RMatrix=[]
    for i in range(len(trs)):
        totalrow=m[trs[i]]
        # print(totalrow)
        Rtemp=[]
        for j in range(len(abso)):
            Rtemp.append(totalrow[abso[j]])
        RMatrix.append(Rtemp)
    return RMatrix

#Qmatrix producer
def getQMatrix(trs,asbo):
    RMatrix=[]
    for i in range(len(trs)):
        totalrow=m[trs[i]]
        # print(totalrow)
        Rtemp=[]
        for j in range(len(trs)):
            Rtemp.append(totalrow[trs[j]])
        RMatrix.append(Rtemp)
    return RMatrix

#get matrix by user
def GetTransitionMAtrix():
    Dim = input("Enter the matrix dimension:\n")
    transisionMatrix=[]
    for i in range(int(Dim)):
        temp=[]
        for j in range(int(Dim)):
            print("Enter the value of row",i,"Col",j)
            val=input("")
            temp.append(int(val))
        transisionMatrix.append(temp)
    return transisionMatrix

def fundamentalMatrix(matrix):
    #N=(i-Q)^-1
    transinetNodes=getTransientNodes(matrix)
    absorbingNodes=getAbosrbingNodes(matrix)
    Qmatrix=getQMatrix(transinetNodes,absorbingNodes)
    identity=identityMatrix(len(Qmatrix))
    # print("len(identity)")
    # print(len(identity))
    # print("len(Qmatrix)")
    # print(len(Qmatrix))
    IminusQ=subtract(identity,Qmatrix)
    N=getMatrixInverse(IminusQ)
    return N

     
NodeA=[0.5,0.5,0,0]
NodeB=[0,0.5,0.5,0]
NodeC=[0.5,0,0,0.5]
NodeD=[0,0,0,1]

m=[NodeA,NodeB,NodeC,NodeD]

Fund=fundamentalMatrix(m)
print("Fundemental Matrix: ")
print(Fund)