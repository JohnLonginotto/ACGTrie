from distutils.core import setup, Extension
from Cython.Build import cythonize


setup(
    packages=['acgtrie'],
    ext_modules=cythonize(Extension(
        'acgtrie/cpp_acgtrie',
        sources=['acgtrie/cpp_acgtrie.pyx', 'src/cpp_acgtrie.cpp'],
        extra_compile_args=['-std=c++11', '-O3', '-march=native'],
        extra_link_args=['-lstdc++'],
        language='c++',
    )),
)
