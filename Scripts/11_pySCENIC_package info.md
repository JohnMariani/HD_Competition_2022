```python
import os
import glob
import pickle 
import pandas
import numpy
import dask
import arboreto
import pyscenic
```


```python
print('\n'.join(f'{m.__name__}=={m.__version__}' for m in globals().values() if getattr(m, '__version__', None)))
```

    pandas==1.2.2
    numpy==1.17.0
    dask==2021.02.0
    pyscenic==0.11.0



```python
print("arboreto==0.1.6")
```

    arboreto==0.1.6

