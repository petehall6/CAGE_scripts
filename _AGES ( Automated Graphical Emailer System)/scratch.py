import sys
import os

functions = os.path.join(os.getcwd(),'functions')

print(functions)

sys.path.append(functions)

import test_script


a = test_script.found_it()

a