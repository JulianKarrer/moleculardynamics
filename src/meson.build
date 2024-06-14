lib_sources = [  # All source files (excluding headers)
    'hello.cpp',
    'verlet.cpp',
    'atoms.cpp',
    'lj_direct_summation.cpp',
    'thermostat.cpp',
    'neighbors.cpp',
    'ljts.cpp',
    'xyz.cpp'
]

lib = static_library(
    'lib',  # Unique name
    lib_sources,
    dependencies : [eigen, mpi]
)

lib_incdirs = include_directories('.')
