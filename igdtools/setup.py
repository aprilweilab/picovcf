import os
import subprocess
import sys
import shutil
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

THIS_DIR = os.path.dirname(os.path.realpath(__file__))

env_debug = int(os.environ.get("IGDTOOLS_DEBUG", 0))
env_no_vcfgz = int(os.environ.get("IGDTOOLS_NO_VCFGZ", 0))

# This is for handling `python setup.py bdist_wheel`, etc.
extra_cmake_args = []
build_type = "Release"
enable_vcfgz = "ON"

if env_debug:
    build_type = "Debug"
if env_no_vcfgz:
    enable_vcfgz = "OFF"

# This is a bit of hack. We want to always have the up to date picovcf, and the source dist
# for igdtools copies the relevent contents of _this_ directory. So any time this file is
# parsed we just copy.
if os.path.isfile("../picovcf.hpp"):
    shutil.copy2("../picovcf.hpp", "picovcf.hpp")

# igdtools is not really a Python C extension, it is just an executable that is wrapped
# with Python code so that it can be pip installed. So the resulting Wheel should not be
# parameterized by Python version, just 
# a Wheel without using the build_ext stuff though, so we have to specialize to Python versions
# even though technically it is not needed from a platform point of view.

class CMakeExtension(Extension):
    def __init__(self, name, cmake_lists_dir=".", sources=[], extra_executables=[], **kwa):
        Extension.__init__(self, name, sources=sources, **kwa)
        self.cmake_lists_dir = os.path.abspath(cmake_lists_dir)
        self.extra_executables = extra_executables

class CMakeBuild(build_ext):
    def get_source_files(self):
        return ["CMakeLists.txt", "picovcf.hpp", "igdtools.cpp", "third-party/args.hxx"]

    def build_extensions(self):
        assert len(self.extensions) == 1
        ext = self.extensions[0]
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        try:
            subprocess.check_call(["cmake", "--version"])
        except OSError:
            raise RuntimeError("Cannot find CMake executable")

        cmake_args = [
            "-DCMAKE_BUILD_TYPE=%s" % build_type,
            "-DCMAKE_RUNTIME_OUTPUT_DIRECTORY_{}={}".format(build_type.upper(), "."),
            "-DENABLE_VCF_GZ=%s" % enable_vcfgz,
        ] + extra_cmake_args
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)


        # Config and build the executable
        subprocess.check_call(["cmake", ext.cmake_lists_dir] + cmake_args,
                                cwd=self.build_temp, stdout=sys.stdout)
        subprocess.check_call(["cmake", "--build", ".", "--config", build_type, "--", "-j"],
                                cwd=self.build_temp, stdout=sys.stdout)
        #shutil.copy(f"{BUILD_DIR}/{BINARY_NAME}", f"{PACKAGE_NAME}/{BINARY_NAME}")
        for executable in ext.extra_executables:
            shutil.copy2(os.path.join(self.build_temp, executable),
                         os.path.join(extdir, executable))

PACKAGE_NAME = "igdtools"
BINARY_NAME = "igdtools"
version = "2.1"
with open(os.path.join(THIS_DIR, "README.md")) as f:
    long_description = f.read()
setup(
    name=PACKAGE_NAME,
    packages=[PACKAGE_NAME],
    version=version,
    description="Tools for converting VCF to IGD files and processing them.",
    author="Drew DeHaas",
    author_email="",
    url="https://github.com/aprilweilab/picovcf",
    zip_safe=False,
    ext_modules=[CMakeExtension(BINARY_NAME, extra_executables=[BINARY_NAME])],
    cmdclass={"build_ext": CMakeBuild},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: MacOS",
        "Operating System :: POSIX :: Linux",
    ],
    entry_points={
        "console_scripts": ["igdtools=igdtools.main:main"],
    },
    long_description=long_description,
    long_description_content_type="text/markdown",
)
