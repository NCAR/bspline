# -*- python -*-

import eol_scons
from SCons.Script import Environment

def bspline_build(env):
    "Global tool for building just bspline."
    env.Require('buildmode')
    env.Append(CXXFLAGS=['-Wextra', '-std=c++11'])

env = Environment(tools=['default'], GLOBAL_TOOLS=[bspline_build])

# The primary builds are in the tool.  This file adds the testing products,
# for when the bspline is source is built by itself.

env.SConscript('BSpline/tool_bspline.py')
env.SConscript('Tests/C++/SConscript')
