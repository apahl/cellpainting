## Compile this Nim module using the following command:
##   python3 path/to/pmgen.py mod_profile_sim.nim
## Then copy the generated `nim_ext.so` (Linux) to the cellpainting folder.

import pymod
import pymodpkg/docstrings

# parameter type bool is not implemented yet
proc profile_sim*(current, reference: string; ignore_direction: int=0; mask: string=""): float {.exportpy.} =
  docstring"""Calculate the similarity of two activity_profiles of the same length.
  ignore_direction encodes a bool as int (bool not yet implemented by pymod), 0 is false.
  When mask is given, only the significant in the mask are taken into account.
  Returns value between 0 .. 1"""
  let
    # ignore_direction = 0
    ignore_dir = if ignore_direction == 0: false else: true
    ref_len = len(reference)
    mask_len = len(mask)
  doAssert(ref_len == len(current), "Activity Profiles must have the same length to be compared.")
  var
    matching_bytes = 0
    sig_pos = 0
  for idx in 0 .. ref_len-1:
    # when a mask is used only consider those positions that are significant in the mask
    if mask_len > 0 and mask[idx] == '1': continue
    if current[idx] != '1' or reference[idx] != '1':
      sig_pos += 1
      if ignore_dir:
        if current[idx] != '1' and reference[idx] != '1':
          matching_bytes += 1
      else:
        if current[idx] == reference[idx]:
          matching_bytes += 1
  result = matching_bytes / sig_pos

initPyModule("nim_ext", profile_sim)
