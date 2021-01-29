# -*- python -*-

import eol_scons
from SCons.Script import Environment, Copy

env = Environment(tools=['default'])

env.Export('env')

env.SConscript('SConscript')
env.SConscript('Tests/C++/SConscript')

IMAGES = env.Split("""
plot-1.png plot-2.png plot-3.png plot-4.png plot-5.png plot-6.png
""")

IMAGES = [env.File("#/Tests/R/" + f) for f in IMAGES]

docfiles = env.Split("""
 README.md
 Doxyfile
 BSpline/BSpline.dox
 BSpline/BSplineLib.cpp
 BSpline/BSpline.cpp
 BSpline/BSpline.h
 BSpline/BSplineBase.cpp
 BSpline/BSplineBase.h
 BSpline/BandedMatrix.h
""")
docfiles.extend([str(f) for f in IMAGES])

doc = env.Command(env.File('#/doc/index.html'), docfiles,
                  ['doxygen'] +
                  [Copy(env.Dir('#/doc'), f) for f in IMAGES])

env.Alias('doc', doc)
