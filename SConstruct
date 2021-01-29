# -*- python -*-

import eol_scons
from SCons.Script import Environment, Copy

env = Environment(tools=['default', 'doxygen'])

env.Export('env')

env.SConscript('SConscript')
env.SConscript('Tests/C++/SConscript')

IMAGES = env.Split("""
plot-1.png plot-2.png plot-3.png plot-4.png plot-5.png plot-6.png
""")

IMAGES = [env.File("#/Tests/R/" + f) for f in IMAGES]

docfiles = env.Split("""
 README.md
 BSpline/BSpline.dox
 BSpline/BSplineLib.cpp
 BSpline/BSpline.cpp
 BSpline/BSpline.h
 BSpline/BSplineBase.cpp
 BSpline/BSplineBase.h
 BSpline/BandedMatrix.h
""")
docfiles.extend([str(f) for f in IMAGES])

doxyfile_text = """
EXTRACT_ALL	       = NO
EXTRACT_STATIC	       = NO
EXTRACT_PRIVATE	       = NO

SOURCE_BROWSER         = NO
REFERENCED_BY_RELATION = NO
CLASS_GRAPH            = YES
COLLABORATION_GRAPH    = NO
INCLUDE_GRAPH          = NO
INCLUDED_BY_GRAPH      = NO
ALPHABETICAL_INDEX     = NO

GRAPHICAL_HIERARCHY    = YES
"""

version = "%(REPO_REVISION)s" % env
docs = env.Apidocs(docfiles,
                   DOXYFILE_TEXT=doxyfile_text,
                   DOXYFILE_DICT={'PROJECT_NAME': 'EOL BSpline Library',
                                  'PROJECT_NUMBER': version})

env.AddPostAction(docs, [Copy(docs[0].dir, f) for f in IMAGES])

env.Alias('doc', docs)
