from distutils.core import setup, Extension

lpn_module = Extension('randLPN_CPP',sources = ['randLPN_CPP.cpp'])
setup(name = 'randLPN_CPP',version='0.0',description = 'Generate random LPNs.',ext_modules = [lpn_module])

