#### Package info for 12_DCA_Denoise.ipynb


```python
import warnings
import numpy
import matplotlib
import scanpy
import pandas
import anndata
import dca
warnings.simplefilter(action="ignore", category = FutureWarning)
```

```python
print('\n'.join(f'{m.__name__}=={m.__version__}' for m in globals().values() if getattr(m, '__version__', None)))
```

    numpy==1.21.5
    matplotlib==3.5.1
    scanpy==1.4.4.post1
    pandas==1.4.0
    anndata==0.6.22.post1



```python
print("matplotlib==3.5.1")
print("dca==0.3.4")
```

    matplotlib==3.5.1
    dca==0.3.4

