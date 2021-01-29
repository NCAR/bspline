
from SCons.Script import Environment, Copy, Export

env = Environment(tools=['default', 'doxygen', 'gitinfo'])

# Since this is often a submodule, make sure to use the git version info for
# this subdirectory.
env.LoadGitInfo('.')
version = env.get('REPO_REVISION', 'unknown')

# If git revision info is loaded successfully, then set the auto-revision,
# otherwise the default will be used.  Someday the default version could be
# extracted from BSplineVersion.h so it can be passed to doxygen.
if version != 'unknown':
    env.Append(CPPDEFINES=['BSPLINE_AUTO_REVISION="%s"' % (version)])

sources = ['BSplineLib.cpp']

bsplinelib = env.Library('bspline', sources)

env.Default(bsplinelib)

tooldir = env.Dir('.')


def bspline(env):
    # The library is static, so add it directly to LIBS.
    env.Append(LIBS=[bsplinelib])
    # Includes are qualified with the BSpline directory, so add the parent
    # to CPPPATH.
    env.AppendUnique(CPPPATH=[tooldir.Dir('..')])


Export('bspline')


# Setup apidocs

IMAGES = env.Split("""
plot-1.png plot-2.png plot-3.png plot-4.png plot-5.png plot-6.png
""")

IMAGES = [env.File("../Tests/R/" + f) for f in IMAGES]

docfiles = env.Split("""
 ../README.md
 BSpline.dox
 BSplineLib.cpp
 BSpline.cpp
 BSpline.h
 BSplineBase.cpp
 BSplineBase.h
 BandedMatrix.h
""")
docfiles.extend(IMAGES)

doxyfile_text = """
EXTRACT_ALL            = NO
EXTRACT_STATIC         = NO
EXTRACT_PRIVATE        = NO

SOURCE_BROWSER         = NO
REFERENCED_BY_RELATION = NO
CLASS_GRAPH            = YES
COLLABORATION_GRAPH    = NO
INCLUDE_GRAPH          = NO
INCLUDED_BY_GRAPH      = NO
ALPHABETICAL_INDEX     = NO

GRAPHICAL_HIERARCHY    = YES
"""

docs = env.Apidocs(docfiles,
                   DOXYFILE_TEXT=doxyfile_text,
                   DOXYFILE_DICT={'PROJECT_NAME': 'EOL BSpline Library',
                                  'PROJECT_NUMBER': version})

env.AddPostAction(docs, [Copy(docs[0].dir, f) for f in IMAGES])

env.Alias('doc', docs)
