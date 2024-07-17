# epistasis_reanalysis

## Creation of conda environment

```
mamba create --name chimera \
anaconda::networkx \
anaconda::numpy \
anaconda::scipy \
anaconda::seaborn \
anaconda::pandas \
anaconda::statsmodels \
conda-forge::matplotlib-venn \
anaconda::ipykernel

mamba activate chimera

pip install hypernetx
```