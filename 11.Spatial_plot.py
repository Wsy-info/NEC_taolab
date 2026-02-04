import stereo as st
import warnings
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from distributed.utils import palette

warnings.filterwarnings('ignore')

### NEC1
### cellbins
data_path = 'NEC1.cellbin.gef'
data = st.io.read_gef(
        file_path=data_path,
        is_sparse=True,
        bin_type='cell_bins',
        )
data.tl.cal_qc()

coords_array = data.cell_borders
cell_names = data.cell_names.astype(str).tolist()

meta = pd.read_csv("../NEC1.csv")
meta = meta.rename(columns={'Unnamed: 0': 'cell_name'})

valid_mask = (coords_array != 32767).all(axis = 2)
valid_coords = coords_array[valid_mask]

cell_indices = np.repeat(np.arange(len(cell_names)), 32)
valid_cell_indices = cell_indices[valid_mask.flatten()]
valid_cell_names = [cell_names[i] for i in valid_cell_indices]

df = pd.DataFrame({
    'cell_name': valid_cell_names,
    'x': valid_coords[:, 0],
    'y': valid_coords[:, 1]
})

df['cell_name'] = 'NEC1_' + df['cell_name']
df_meta = df.merge(meta, on='cell_name', how='left')
df_meta.to_csv("../NEC1_all.csv")