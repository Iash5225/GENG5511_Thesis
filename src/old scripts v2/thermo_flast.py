import ctypes
dll_path = r"C:\Users\iashb\AppData\Roaming\Mathematica\Paclets\Repository\ThermoFAST64-2.2.0.1\Windows\ThermoFAST64.dll"

try:
    tf = ctypes.CDLL(dll_path)
    print("DLL loaded successfully!")

    # Try printing the available symbols (if decorated with __declspec(dllexport))
    # This will likely error if we don’t know the function name
    print(dir(tf))

except Exception as e:
    print("Error loading DLL:", e)
