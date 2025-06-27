import platform
import subprocess
from setuptools import setup, Extension

# Function to check for Fortran compiler


def has_fortran_compiler():
    try:
        # Try a common Fortran compiler command
        subprocess.check_output(['gfortran', '--version'], stderr=subprocess.STDOUT)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        try:
            subprocess.check_output(['ifort', '--version'], stderr=subprocess.STDOUT)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            return False


extensions = []

if has_fortran_compiler():
    flags = []
    libs = []
    if platform.system() == 'Windows':
        flags = []
        libs = []
    elif platform.system() == 'Darwin':
        flags = ['-fPIC']
        libs = ['/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib']
    else:
        flags = ['-fPIC']
        libs = []

    define_macros = [('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')]

    extensions = [
        Extension('seapy.external.oalib', sources=['src/oalib.f'],
                  extra_f77_compile_args=flags,
                  library_dirs=libs,
                  define_macros=define_macros),
        Extension('seapy.external.hindices', sources=['src/hindices.f'],
                  extra_f77_compile_args=flags,
                  library_dirs=libs,
                  define_macros=define_macros),
        Extension('seapy.external.extractobs', sources=['src/extract_obs.f'],
                  extra_f77_compile_args=flags,
                  library_dirs=libs,
                  define_macros=define_macros),
    ]
else:
    print("WARNING: No FORTRAN compiler found. Seapy Extensions will not be built.")

setup(
    ext_modules=extensions,
    # No other metadata here, it's all in pyproject.toml
)
