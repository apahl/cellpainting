# Nim Extensions for Speedup

This folder contains faster versions of some Python functions used in this project, written in [Nim](https://nim-lang.org/). They can be compiled into Python extensions using the excellent [Pymod](https://github.com/jboy/nim-pymod) module.
The use of these extensions is completely optional, the project uses the pure Python versions when it can not find the extensions.


Compile this Nim module with the following command (requires a working Nim toolchain):
`python3 path/to/pmgen.py mod_profile_sim.nim`.
Then copy the generated `nim_ext.so` (Linux, not tested on Windows) to the cellpainting folder.