from distutils.extension import Extension

ext_modules = [Extension("mpcincfgm",
    ["mpcinc/cython/mpcincfgm.pyx",
        "mpcinc/mpcincfgmdynmem.c","mpcinc/mpcincfgm.c",
        "mpcinc/mpcinccvpdynmem.c","mpcinc/mpcinccvp.c",
        "mpcinc/cjson.c", "mpcinc/mpcincdynmem.c",
        "mpcinc/mpcincmtxops.c"],
    include_dirs = ["mpcinc/include"], library_dirs = ["mpcinc"],
    libraries=["m"])]
