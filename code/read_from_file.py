""" Read CSI information from file, logged by log_to_file.c tool

Astrit Zhushi (c) 2011 a.zhushi@cs.ucl.ac.uk

"""
from iwlnl_struct import *


def read_from_file(csi_log_file):

    with open(csi_log_file, 'rb') as f:
        data = f.read()

    cur = 0
    length = len(data)
    csi = []

    while cur < (length - 3):
        (field_length,) = struct.unpack(">H", data[cur:cur + 2])
        (code,) = struct.unpack("B", data[cur + 2:cur + 3])
        cur += 3

        csi_bytes = data[cur:cur + field_length - 1]
        cur = cur + field_length - 1

        if code != 187:
            print("Unhandled code %d " % code)
            continue
        iwlnstruct = iwlnl_struct(csi_bytes, from_file=True)
        csi.append(iwlnstruct)
    return csi
