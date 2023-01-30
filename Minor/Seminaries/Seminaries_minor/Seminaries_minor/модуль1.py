import functools



class Tensor:
    def __init__(self, data, shape):
        self.data = data
        self.shape = shape

    def __add__(self, other):
        if self.shape != other.shape:
            raise ValueError("Tensors must have the same shape for addition.")
        return Tensor([a + b for a, b in zip(self.data, other.data)], self.shape)

    def __mul__(self, other):
        if isinstance(other, (int,float)):
            return Tensor([x * other for x in self.data], self.shape)
        elif isinstance(other, Tensor):
            if self.shape[-1] != other.shape[0]:
                raise ValueError("Incompatible tensors for contraction")
            new_shape = self.shape[:-1] + other.shape[1:]
            new_data = sum(elem * other_elem for elem, other_elem in zip(self.elements, other.elements))
            return Tensor(new_data,new_shape)
        else:
            raise ValueError("Wrong data type")
            
    def __str__(self):
        return f"Tensor with shape {self.shape} and data {self.data}"

    def contract(self, other, *indices):
        if len(indices) != 2:
            raise ValueError("Exactly two indices must be contracted.")
        if indices[0] >= len(self.shape) or indices[1] >= len(other.shape):
            raise ValueError("Invalid indices for contraction.")
        if self.shape[indices[0]] != other.shape[indices[1]]:
            raise ValueError("Tensors must have matching dimensions at contracted indices.")

        new_shape = self.shape[:indices[0]] + self.shape[indices[0]+1:] + other.shape[:indices[1]] + other.shape[indices[1]+1:]
        new_data = []
        for i in range(len(self.data)):
            for j in range(len(other.data)):
                si = i // functools.reduce(lambda x, y: x*y, self.shape[indices[0]+1:])
                sj = j // functools.reduce(lambda x, y: x*y, other.shape[indices[1]+1:])
                if si % self.shape[indices[0]] == sj % other.shape[indices[1]]:
                    new_data.append(self.data[i] * other.data[j])
        return Tensor(new_data, new_shape)

# Define two tensors
A = Tensor([[[1, 2], [3, 4]], [[5, 6], [7, 8]]], shape=(2, 2, 2))
B = Tensor([[[9, 10], [11, 12]], [[13, 14], [15, 16]]], shape=(2, 2, 2))

# Perform tensor multiplication
C = A * B

# Print the resulting tensor
print(C)