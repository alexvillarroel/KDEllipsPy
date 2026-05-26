from setuptools import setup, find_packages
import subprocess
import os
from setuptools.command.install import install
from setuptools.command.develop import develop

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
            subprocess.check_call(['make', 'python'], cwd=axitra_src)
            print("Axitra build process completed successfully.")
        except Exception as e:
            print(f"Error during Axitra build: {e}")
            print("Requirements: gfortran, make, and numpy<2.0.0 must be available in the environment.")
    else:
        print(f"Axitra source directory not found at: {axitra_src}")

class PostInstallCommand(install):
    def run(self):
        install.run(self)
        build_axitra()

class PostDevelopCommand(develop):
    def run(self):
        develop.run(self)
        build_axitra()

setup(
    name="kdellipspy",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.26,<2.0.0",
        "scipy",
        "matplotlib",
        "pyproj",
        "obspy",
        "PySide6",
        "neighpy",
        "pymc>=5.9.0",
        "arviz",
    ],
    entry_points={
        "console_scripts": [
            "kdellipspy=kdellipspy.cli:main",
        ],
    },
    cmdclass={
        'install': PostInstallCommand,
        'develop': PostDevelopCommand,
    },
)
