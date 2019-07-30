import os
import sys

is_windows = os.name == "nt"
is_apple = sys.platform == "darwin"
if is_windows:
    suffix= ".dll"
elif is_apple:
    suffix = ".dylib"
else:
    suffix = ".so"
