# -*- python -*-

import eol_scons
from SCons.Script import Environment

env = Environment(tools=['default'])

# The primary builds are in the tool.  This file adds the testing products,
# for when the bspline is source is built by itself.

env.SConscript('BSpline/tool_bspline.py')
env.SConscript('Tests/C++/SConscript')
