## Compile this Nim module using the following command:
##   python3 path/to/pmgen.py mod_profile_sim.nim
## Then copy the generated `nim_ext.so` (Linux) to the cellpainting folder.

import pymod
import pymodpkg/docstrings

proc profile_sim*(current, reference: string): float {.exportpy.} =
  docstring"""Calculate the similarity of two byte activity_profiles of the same length.
  Returns value between 0 .. 1"""
  let ref_len = len(reference)
  doAssert(ref_len == len(current), "Activity Profiles must have the same length to be compared.")
  var matching_bytes = 0
  for idx in 0 .. ref_len-1:
    if current[idx] == reference[idx]:
      matching_bytes += 1
  result = matching_bytes / ref_len

initPyModule("nim_ext", profile_sim)