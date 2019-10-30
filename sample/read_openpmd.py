
import openpmd_api as api

# example: data handling
import numpy as np

file_name = "./fields.h5"
series = api.Series( file_name, api.Access_Type.read_only)

print(list(series.iterations))

from pprint import pprint
#pprint(vars(series))
#pprint(vars(series.iterations))

i = series.iterations[1];

print("openPMD version: ",
      series.openPMD)

# record
cB = i.meshes["B"]

# record components
cbx = cB["x"]

x_data = cbx.load_chunk()

series.flush()

extent = cbx.shape

print(
    "First values in E_x "
    "of shape: ",
    extent)


print(x_data)
