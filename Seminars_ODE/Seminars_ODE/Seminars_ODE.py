import numpy as np
import pandas as pd

names = ("id", "n", "Na", "Mg", "Al", "Si", "K", "Ca", "Ba", "Fe", "Type")
data = pd.read_csv("data.csv", header = None, names =  names, index_col = 0)

#print(data.head()) 

print(data[["n", "Na", "Mg"]].head(10))


