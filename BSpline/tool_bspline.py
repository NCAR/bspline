toolname = 'bspline'

thisDir = Dir('.').abspath
upDir   = Dir('./../').abspath

# Define tool. This must match the value of toolname.
def bspline(env):
    env.Append(LIBS = [toolname])
    env.AppendUnique(CPPPATH = [upDir])
    env.AppendUnique(CPPPATH = [thisDir])
    env.AppendUnique(LIBPATH = [thisDir])
    env.Require(['prefixoptions'])
    
Export(toolname)

# Build library
sources = Split("""
BSplineLib.cpp
""")

headers = Split("""
BSplineBase.h
BSpline.h
""")

env = Environment(tools = ['default', 'doxygen', toolname])

lib = env.Library(toolname, sources)
env.Default(lib)

# Create doxygen
doxref = env.Apidocs(sources + headers, DOXYFILE_DICT={'PROJECT_NAME':toolname, 'PROJECT_NUMBER':'1.0'})
