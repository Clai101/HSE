import functools

class array:
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
        if isinstance(a, (int,float)):
            return  [[j * a for j in i] for i in self.arr]
        elif isinstance(a, array):
            if self.size[1] != a.size[0]:
                raise ValueError("Incompatible matrix size for multiplication")
            else:
                return array([[functools.reduce(lambda x, y: x + y, [self.arr[i][k] * a.arr[k][j] for k in range(self.size[1])]) 
                             for j in range(a.size[1])] for i in range(self.size[0])])
        else:
            raise ValueError("Wong data type")

    def __add__(self, a):
        if self.si
        if isinstance(a, (int,float)):
            return [[j + a for j in i] for i in self.arr]
        elif isinstance(a, array):
            if self.size != a.size:
                raise ValueError("Incompatible matrix size for multiplication")
            else:
                return array([[a.arr[i][j] + self.arr[i][j] for j in range(self.size[1])] for i in range(self.size[0])])
        else:
            raise ValueError("Wong data type")

    def __str__(self):
        return "\n".join(["\t".join(map(str,row)) for row in self.arr])

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

print(array([1, 2, 3]) *  3)

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