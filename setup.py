from setuptools import setup, find_packages
import subprocess
import os
from setuptools.command.install import install
from setuptools.command.develop import develop

from setuptools.command.build_py import build_py
from setuptools.command.build_ext import build_ext

def build_axitra():
    print("Building axitra binaries and python wrappers...")
    root_dir = os.path.abspath(os.path.dirname(__file__))
    axitra_src = os.path.join(root_dir, 'kdellipspy', 'axitra', 'src')
    
    if os.path.exists(axitra_src):
        try:
            import sys
            # Clean previous builds
            subprocess.check_call(['make', 'clean'], cwd=axitra_src)
            # Build core binaries (axitra, convms)
            print("Compiling Fortran binaries (axitra, convms)...")
            subprocess.check_call(['make', 'all'], cwd=axitra_src)
            # Build python wrappers and utilities (convmPy, axitra2py)
            print("Building Python wrappers (f2py) and utilities...")
            subprocess.check_call(['make', 'python', f'PYTHON={sys.executable}'], cwd=axitra_src)
            print("Axitra build process completed successfully.")
        except Exception as e:
            print(f"Error during Axitra build: {e}")
            print("Requirements: gfortran, make, and numpy>=2.0.0 must be available in the environment.")
    else:
        print(f"Axitra source directory not found at: {axitra_src}")

class CustomBuildPyCommand(build_py):
    def run(self):
        build_axitra()
        super().run()

class CustomBuildExtCommand(build_ext):
    def run(self):
        build_axitra()
        super().run()

setup(
    name="kdellipspy",
    version="0.1.0",
    python_requires=">=3.10, <3.14",
    packages=find_packages(),
    install_requires=[
        "numpy>=2.0.0",
        "scipy>=1.13.0",
        "matplotlib",
        "cartopy",
        "pyproj",
        "obspy",
        "neighpy",
        "joblib",
        "pymc>=5.10.0",
        "arviz>=0.18.0",
    ],
    entry_points={
        "console_scripts": [
            "kdellipspy=kdellipspy.cli:main",
            "kde-prep=kdellipspy.prep_data:main",
            "kde-traces=kdellipspy.plot_traces:main",
            "kde-map=kdellipspy.map:main",
            "kde-dashboard=kdellipspy.dashboard:main",
        ],
    },
    cmdclass={
        'build_py': CustomBuildPyCommand,
        'build_ext': CustomBuildExtCommand,
    },
)
