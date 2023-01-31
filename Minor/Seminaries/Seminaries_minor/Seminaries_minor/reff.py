import functools

class array:
<<<<<<< HEAD
    def __init__(self, arr):
        self.arr = arr
        self.size = []
        while type(arr) == type([0]):
            self.size.append(len(arr))
            arr = arr[0]
        if len(self.size) == 1:
            self.size = [1] + self.size
        for i in range(0,len(self.arr)):
            if self.size[0] != 1:
                if len(self.arr[i]) != self.size[1]:
                    raise ValueError("Input matrix is not a rectangle")
            elif isinstance(self.arr[i], list):
                self.arr[i] = array(self.arr[i])

    def __mul__(self, a):
        if self.size[0] == 1:
            if isinstance(a, (int,float)):
                return array([i * a for i in self.arr])
            elif isinstance(a, array):
                if self.size[1] != a.size[0]:
                    raise ValueError("Incompatible matrix size for multiplication")
                elif self.size[0] == a.size[1] == 1:
                    return float([[functools.reduce(lambda x, y: x + y, [self.arr[i] * a.arr[k][j] for k in range(self.size[1])]) 
                                 for j in range(a.size[1])] for i in range(self.size[0])][0][0])
                else:
                    return array([[functools.reduce(lambda x, y: x + y, [self.arr[i] * a.arr[k][j] for k in range(self.size[1])]) 
                                 for j in range(a.size[1])] for i in range(self.size[0])])
            else:
                raise ValueError("Wong data type")
        else:
            if isinstance(a, (int,float)):
                return  array([[j * a for j in i] for i in self.arr])
            elif isinstance(a, array):
                if self.size[1] != a.size[0]:
                    raise ValueError("Incompatible matrix size for multiplication")
                else:
                    return array([[functools.reduce(lambda x, y: x + y, [self.arr[i][k] * a.arr[k][j] for k in range(self.size[1])]) 
                                 for j in range(a.size[1])] for i in range(self.size[0])])
            else:
                raise ValueError("Wong data type")

    def __add__(self, a):
        if self.size[0] == 1:
            if isinstance(a, (int,float)):
                return array([i + a for i in self.arr])
            elif isinstance(a, array):
                if self.size != a.size:
                    raise ValueError("Incompatible matrix size for multiplication")
                else:
                    return array([a.arr[i] + self.arr[i] for i in range(self.size[0])])
            else:
                raise ValueError("Wong data type")
        else:
            if isinstance(a, (int,float)):
                return array([[j + a for j in i] for i in self.arr])
            elif isinstance(a, array):
                if self.size != a.size:
                    raise ValueError("Incompatible matrix size for multiplication")
                else:
                    return array([[a.arr[i][j] + self.arr[i][j] for j in range(self.size[1])] for i in range(self.size[0])])
            else:
                raise ValueError("Wong data type")
    def __radd__(self, a):
        return array(self.arr) + a
    def __rmul__(self, a):
        return array(self.arr) * a

    def __str__(self):
        if self.size != [1, 1]:
            return "\n".join(["\t".join(map(str,row)) for row in self.arr])
        return str(self.arr[0])
    def Gaussian_elimination(self):
        A = self.arr
        for i in range(len(A)):
            pivot = A[i][i]
            A[i] = A[i] * pivot**(-1)
            for j in range(i+1, len(A)):
                factor = A[j][i]
                A[j] = A[j] - factor * A[i]
        return A



A = array([[2, 5, 5], [2, 0, 2], [1, -0.5, 1]])

print( array([1, 2, 3]) * A)

"""
# Coefficient matrix
A = [[3, 2, -1], [2, -2, 4], [-1, 0.5, -1]]

# Constant term
b = [1, -2, 0]

# Perform Gaussian elimination
for i in range(len(A)):
    # Divide the row by the pivot element
    pivot = A[i, i]
    A[i] = A[i] / pivot
    b[i] = b[i] / pivot
    for j in range(i+1, len(A)):
        # Subtract the row with the pivot element from the other rows
        factor = A[j, i]
        A[j] = A[j] - factor * A[i]
        b[j] = b[j] - factor * b[i]

# Perform back-substitution
x = np.zeros(len(A))
x[-1] = b[-1]
for i in range(len(A)-2, -1, -1):
    x[i] = b[i] - np.dot(A[i], x)[i]

print(x)
"""
=======
    def __init__(self, arr, t = 1):
        if isinstance(arr, (int, float)):
            self.arr, self.size = arr, [1, 1]
        elif isinstance(arr[0], (int, float)) and t < 2:
            t += 1
            self.arr, self.size = [arr], []
            arr = [arr]
            while not isinstance(arr, (int, float)):
                self.size.append(len(arr))
                arr = arr[0]
        else:
            t += 1
            self.arr, self.size = arr, []
            while not isinstance(arr, (int, float)):
                self.size.append(len(arr))
                arr = arr[0]

    def __getitem__(self, indices):
        if isinstance(indices, int):
            return self.arr[indices]
        elif isinstance(indices, tuple):
            if isinstance(indices[0], int) and isinstance(indices[1], int):
                return self.arr[indices[0]][indices[1]]
            elif isinstance(indices[0], slice) and isinstance(indices[1], int):
                return array([row[indices[1]] for row in self.arr[indices[0]]])
            elif isinstance(indices[0], int) and isinstance(indices[1], slice):
                return array([self.arr[indices[0]][i] for i in range(len(self.arr[indices[0]]))[indices[1]]])
            elif isinstance(indices[0], slice) and isinstance(indices[1], slice):
                return array([[self.arr[i][j] for j in range(len(self.arr[i]))[indices[1]]] 
                              for i in range(len(self.arr))[indices[0]]])
        else:
            raise ValueError("Invalid index")

    def __mul__(self, a):
        if isinstance(a, (int,float)):
            return  array([[j * a for j in i] for i in self.arr])
        elif isinstance(a, array):
            if self.size[1] != a.size[0]:
                raise ValueError("Incompatible matrix size for multiplication")
            else:
                return array([[functools.reduce(lambda x, y: x + y, [self.arr[i][k] * a.arr[k][j] for k in range(self.size[1])]) 
                             for j in range(a.size[1])] for i in range(self.size[0])])
        else:
            raise ValueError("Wong data type")

    def __add__(self, a):
        if isinstance(a, (int,float)):
            return array([[j + a for j in i] for i in self.arr])
        elif isinstance(a, array):
            if self.size != a.size:
                raise ValueError("Incompatible matrix size for multiplication")
            else:
                return array([[a.arr[i][j] + self.arr[i][j] for j in range(self.size[1])] for i in range(self.size[0])])
        else:
            raise ValueError("Wong data type")

    def __radd__ (self, a):
        return array(self.arr) + a

    def __rmul__ (self, a):
        return array(self.arr) * a

    def __str__(self):
        return "\n".join(["\t".join(map(str,row)) for row in self.arr])

    @property
    def Gaussian_elimination(self):
        A = self * 1
        n, m = A.size[0], A.size[1]
        for i in range(n):
            if A[i, min(i,m-1)] == 0:
                continue
            A.arr[i] = (A[i:i+1, 0:m]*A[i, min(i, m-1)]**-1).arr[0]
            for j in range(i+1, n):
                A.arr[j] = (A[j:j+1, 0:m] + (-1)*A[i:i+1, 0:m]*A[j, min(i,m-1)]).arr[0]
        for i in range(n-1, 0, -1):
            if A[i, min(i, m-1)] == 0:
                continue
            for j in range(0, i):
                A.arr[j] = (A[j:j+1, 0:m] + (-1)*A[i:i+1, 0:m]*A[j, min(i,m-1)]).arr[0]
        return A
    
    @property
    def determinant(self):
        matrix = self.arr
        if len(matrix) == 2:
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
        else:
            result = 0
            for c in range(len(matrix)):
                minor = array([row[:c] + row[c+1:] for row in matrix[1:]])
                result += ((-1) ** c) * matrix[0][c] * minor.determinant
            return result

    @property
    def decomposition(self):
        A = self * 1
        Ag = A.Gaussian_elimination
        a = []
        for i in range(Ag.size[1]):
            _ = -1
            a = a + [_]
            for j in range(Ag.size[0]):
                if Ag[j, i] == 1 and j > max(a):
                    a[-1] = j
        return array([Ag.arr[i] for i in a if i >= 0]), array([[A.arr[i][j] for j in range(A.size[1]) if a[j] >= 0] for i in range(A.size[0])])
    


A1 = array([[1, 2, 3], [4, 5, 6], [1, 2, 4]])
A2 = array([[2, -1, 0], [-1, 1, 1], [0, 1, 2]])
A3 = array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
A4 = array([[1, 0, 0], [1, 1, 1], [0, 1, 1]])
A5 = array([[1, 0, 0], [1, 0, 1], [0, 0, 1]])
A6 = array([[4, 2, 8, 5], [7, 45, 2, 0], [76, 1, 3, 3], [5, 6, 4, 7]])

P1 = array([[1, 0], [1, 0], [0, 3]])

A = P1


F, G = A.decomposition
print(F, G, sep = '\n\n\n')
<<<<<<< HEAD
print(G * F, sep = '\n\n\n')
=======
print(G * F, sep = '\n\n\n')
>>>>>>> 11d9c2dffe39db0a3f04a1c838ba68a75f4ae955
>>>>>>> b7e132ce7b05900c01b0f055621f435bf49336f4
