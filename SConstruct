#
# SCons build script for midapack/mappraiser
#

# import SCons
import os
import subprocess

# import platform
# import re

############################################################
#
# Compiler, compilation flags, etc.
#

cc = "cc"

cflags = Split("-fPIC -xcore-avx2 -axmic-avx512 -pthread -std=gnu99")
optim = Split("-O3")
debug = Split("-g")
# openmp   = Split( '-qopenmp' )
defines = Split("W_MPI")
wecg = Split("W_ECG")
wopenmp  = Split("W_OPENMP")

############################################################
#
# define environment
#

opts = Variables(None, ARGUMENTS)

opts.Add(
    PathVariable(
        "prefix", "installation directory", ".", validator=PathVariable.PathIsDirCreate
    )
)

opts.Add(BoolVariable("fullmsg", "set to enable full output", 0))
opts.Add(BoolVariable("debug", "set to enable building with debug informations", 0))
opts.Add(BoolVariable("optimise", "set to enable building with optimisation", 1))
# opts.Add( BoolVariable( 'profile',  'set to enable building with profile informations', 0 ) )

opts.Add(BoolVariable("shared", "set to enable building of shared libraries", 1))
opts.Add(BoolVariable("static", "set to enable building of static libraries", 0))

opts.Add(BoolVariable("openmp", "set to enable building with OpenMP activated", 0))
opts.Add(BoolVariable("ecg", "set to enable ECG solver", 0))

env = Environment(
    options=opts,
    ENV=os.environ,
    CC=cc,
    CFLAGS=cflags,
    CPPDEFINES=defines,
    #    LIBS       = ['m']
)

if not env["fullmsg"]:
    env.Replace( CCCOMSTR        = " CC     $SOURCES" )
    env.Replace( SHCCCOMSTR      = " CC     $SOURCES" )
    env.Replace( CXXCOMSTR       = " C++    $SOURCES" )
    env.Replace( SHCXXCOMSTR     = " C++    $SOURCES" )
    env.Replace( FORTRANCOMSTR   = " FC     $SOURCES" )
    env.Replace( SHFORTRANCOMSTR = " FC     $SOURCES" )
    env.Replace( LINKCOMSTR      = " Link   $TARGET" )
    env.Replace( SHLINKCOMSTR    = " Link   $TARGET" )
    env.Replace( ARCOMSTR        = " AR     $TARGET" )
    env.Replace( RANLIBCOMSTR    = " Index  $TARGET" )
    env.Replace( INSTALLSTR      = " Install   $TARGET" )

if env["optimise"]:
    env.Append(CFLAGS=optim)

if env["debug"]:
    env.Append(CFLAGS=debug)
    
if env["openmp"]:
    env.Append(CPPDEFINES=wopenmp)

# static overrides shared
if env["static"]:
    env["shared"] = 0

if os.environ.get("NERSC_HOST", None) is not None:
    # we are at NERSC, use Cray compiler wrappers
    cray_opts = subprocess.check_output(["cc", "--cray-print-opts"])
    d = env.ParseFlags(cray_opts)
    env.MergeFlags(d)

if env.GetOption("clean"):
    print("Skip environment conf checks")
else:
    conf = Configure(env)

    # Checks for headers, libraries, etc.
    if not conf.CheckCC():
        Exit(1)
    if not conf.CheckLibWithHeader("mpich", "mpi.h", language="c"):
        # print('Did not find MPI, exiting!')
        Exit(1)
    if not conf.CheckLibWithHeader("fftw3", "fftw3.h", language="c"):
        # print('Did not find FFTW3, exiting!')
        Exit(1)

    env = conf.Finish()

# print("CC is: %s" % env["CC"])
# print("CFLAGS is: %s" % env["CFLAGS"])

# create environments for modules
mapmat_env = env.Clone()
toeplitz_env = env.Clone()
mappraiser_env = env.Clone()

if mappraiser_env['ecg']:
    mappraiser_env.Append(CPPDEFINES=wecg)

# export environments
Export("env")
Export("mapmat_env")
Export("toeplitz_env")
Export("mappraiser_env")

############################################################
#
# Help Text
#

Help(
"""
Type
  'scons'            to build Midapack.
  'scons options'    to get help on command line build parameters.
  'scons install'    to install Midapack in the chosen directory.
  'scons doc'        to create documentation for Midapack.
  'scons dist'       to create tar archive of Midapack.
"""
)

############################################################
#
# build
#

SConscript("src/SConscript")
