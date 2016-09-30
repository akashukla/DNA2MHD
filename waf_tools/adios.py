"""
Detect ADIOS includes and libraries.
"""

import os, glob, types
from waflib.Configure import conf

def options(opt):
    opt.add_option('--enable-adios', help=('Enable parallel build'),
                   dest='enable_adios', action='store_true',
                   default=False)
    opt.add_option('--adios-inc-dir', type='string', help='Path to ADIOS includes', dest='adiosIncDir')
    opt.add_option('--adios-lib-dir', type='string', help='Path to ADIOS libraries', dest='adiosLibDir')

@conf
def check_adios(conf):
    opt = conf.options
    conf.env['ADIOS_FOUND'] = False

    if not conf.options.enable_adios:
	return
    if conf.options.adiosIncDir:
	conf.env.INCLUDES_ADIOS = conf.options.adiosIncDir
    else:
        conf.fatal("Please specify ADIOS include directories by --adios-inc-dir")
    if conf.options.adiosLibDir:
	conf.env.LIBPATH_ADIOS = conf.options.adiosLibDir
        if conf.options.enable_mpi:
            conf.env.LIB_ADIOS = ["adios"]
        else:
            conf.env.append_value('CXXFLAGS', '-D_NOMPI')
            conf.env.append_value('CFLAGS', '-D_NOMPI')
            conf.env.LIB_ADIOS = ["adios_nompi"]            
    else:
        conf.fatal("Please specify ADIOS library directories by --adios-lib-dir")
        
    conf.start_msg('Checking for ADIOS')
    conf.check(header_name='adios.h', features='cxx cxxprogram', use="ADIOS MPI", mandatory=True)
    conf.end_msg("Found ADIOS")
    conf.env['ADIOS_FOUND'] = True
    return 1

def detect(conf):
    return detect_adios(conf)