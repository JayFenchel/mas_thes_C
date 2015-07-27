from distutils.core import setup
from distutils.extension import Extension

ext_modules = [Extension("mpcinccvp",
    ["mpcinc/cython/mpcinccvp.pyx", "mpcinc/cython/mpcincCcvp.pyx",
        "mpcinc/mpcinccvp.c",
    "mpcinc/mpcincdynmem.c",
    "mpcinc/mpcinccvpdynmem.c",
    "mpcinc/cjson.c",  "mpcinc/mpcincmtxops.c"],
    include_dirs = ["mpcinc/include"], library_dirs = ["mpcinc"],
    libraries=["m"]),
    Extension("mpcincCcvp",
    ["mpcinc/cython/mpcincCcvp.pyx",
    ],
    include_dirs = ["mpcinc/include"], library_dirs = ["mpcinc"],
    libraries=["m"])]

