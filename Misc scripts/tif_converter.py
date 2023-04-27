import os
import pandas as pd
from PIL import Image

dir=os.getcwd()
file = 'gel.tif'
tif = os.path.join(dir,file)
output="gel.jpg"
im = Image.open(tif)
im.save(output, "JPEG", quality=100)

