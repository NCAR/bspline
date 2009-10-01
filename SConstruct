# -*- python -*-

import os
import SCons

env = Environment(tools = ['default', 'packaging'])

Export('env')

SConscript('SConscript')
SConscript('Tests/C++/SConscript')

VERSION = str(ARGUMENTS.get('VERSION', ''))
env['PUBDIR'] = "/net/www/docs/homes/granger/bspline"

IMAGES = Split("""
plot-1.png plot-2.png plot-3.png plot-4.png plot-5.png plot-6.png
""")

IMAGES = [ os.path.join("Tests", "R", f) for f in IMAGES ]

distfiles = Split("""
 README
 COPYRIGHT 
 Doxyfile
 SConstruct
 SConscript
 BSpline.sln
 BSpline/BSpline.dox
 BSpline/BSplineLib.cpp
 BSpline/BSpline.cpp
 BSpline/BSpline.h 
 BSpline/BSplineBase.cpp
 BSpline/BSplineBase.h 
 BSpline/BandedMatrix.h
 BSpline/BSpline.vcproj
 BSpline/SConscript
 Tests/Data/sample.temps
 Tests/Data/sample.txt
 Tests/Data/sample.wdir
 Tests/Data/sample.wspd
 Tests/C++/SConscript
 Tests/C++/bspline.cpp
 Tests/C++/options.cpp
 Tests/C++/options.h
""")
distfiles.extend( [ str(f) for f in IMAGES ] )

doc = env.Command(env.File('#/doc/index.html'), distfiles,
                  ['doxygen'] +
                  [ Copy(env.Dir('#/doc'), f) for f in IMAGES ])

env.Alias('doc', doc)

# We have to kludge our own version of this function from the
# SCons.Tools.packaging tool, so that we can properly handle the 'doc'
# directory.  Otherwise scons tries to treat it as a File node and breaks,
# even though the CopyAs action and the Zip tool accept directory nodes.
# Then rather than invoke the Package() builder, invoke the Zip() builder
# directly.
def copytopackageroot(source, env, pkgroot):
    if SCons.Util.is_String(pkgroot):  pkgroot=env.Dir(pkgroot)
    if not SCons.Util.is_List(source): source=[source]

    new_source = []
    for path in source:
        if path.endswith('/'):
            node = env.Dir(path)
            new_node = env.CopyAs(pkgroot.Dir(path), node)
            env.AddPreAction(new_node, Delete(pkgroot.Dir(path)))
        else:
            node = env.File(path)
            new_node = env.CopyAs(pkgroot.File(path), node)
        new_source.extend(new_node)
    return new_source

# If a VERSION has been specified, then create the distribution targets.
if VERSION:

    zipfileroot = "eol-bspline-%s" % VERSION
    zip = env.Zip(zipfileroot, copytopackageroot(distfiles + ['doc/'],
                                                 env,
                                                 zipfileroot))
    env.Alias('zip', zip)

    pubzip = env.Install("$PUBDIR", zip)
    env.AddPostAction(pubzip, "cd $PUBDIR && rm -rf doc")
    env.AddPostAction(pubzip, Copy("$PUBDIR/doc", env.Dir("doc")))
#    env.AlwaysBuild(pubdoc)
#    env.Depends(pubdoc, pubzip)
    env.Alias('publish', pubzip)
#    env.Alias('publish', pubdoc)
    
