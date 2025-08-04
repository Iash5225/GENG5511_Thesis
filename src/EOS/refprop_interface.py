# scripts/refprop_interface.py

import os
import sys
import numpy as np
import comtypes.client


class REFPROP:
    def __init__(self, dll_path='C:\\Refprop\\REFPROP'):
        self.RP = comtypes.client.CreateObject('Refprop.Application')
        self.RP.SETPATHdll(dll_path)
        self.RP.FLAGSdll("Reset HMX", 1)
        self.RP.FLAGSdll("Peng-Robinson", 0)

    def get_triple_point(self, fluid='argon'):
        self.RP.SETFLUIDSdll(fluid)
        T_triple = self.RP.REFPROPdll(
            fluid, "TP", "TTRP", 2, 0, 1, 300, 1, [1.0]).Output[0]
        P_triple = self.RP.REFPROPdll(
            fluid, "TP", "PTRP", 2, 0, 1, 300, 1, [1.0]).Output[0]
        return T_triple, P_triple
