# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
import pyworkflow.em
from pyworkflow.utils import Environ, join

_logo = "nysbc_logo.png"
_references = ['tan2017']

NYSBC_3DFSC_HOME_VAR = 'NYSBC_3DFSC_HOME'


# The following class is required for Scipion to detect this Python module
# as a Scipion Plugin. It needs to specify the PluginMeta __metaclass__
# Some function related to the underlying package binaries need to be
# implemented
class Plugin:
    #__metaclass__ = pyworkflow.em.PluginMeta

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch 3DFSC. """
        environ = Environ(os.environ)
        NYSBC_3DFSC_HOME = os.environ[('%s' % NYSBC_3DFSC_HOME_VAR)]

        environ.update({
            'PATH': join(NYSBC_3DFSC_HOME, 'ThreeDFSC'),
        }, position=Environ.BEGIN)

        if 'PYTHONPATH' in environ:
            # this is required for python virtual env to work
            environ.set('PYTHONPATH', '', position=Environ.BEGIN)
        return environ

    @classmethod
    def getVersion(cls):
        path = os.environ[NYSBC_3DFSC_HOME_VAR]
        for v in cls.getSupportedVersions():
            if v in path:
                return v
        return ''

    @classmethod
    def getSupportedVersions(cls):
        """ Return the list of supported binary versions. """
        return ['2.5']

    @classmethod
    def validateInstallation(cls):
        """ This function will be used to check if package is properly installed. """
        environ = cls.getEnviron()
        missingPaths = ["%s: %s" % (var, environ[var])
                        for var in [NYSBC_3DFSC_HOME_VAR]
                        if not os.path.exists(environ[var])]

        return (["Missing variables:"] + missingPaths) if missingPaths else []

    @classmethod
    def getProgram(cls):
        """ Return the program binary that will be used. """
        if NYSBC_3DFSC_HOME_VAR not in os.environ:
            return None
        cmd = join(os.environ['NYSBC_3DFSC_HOME'], 'ThreeDFSC',
                   'ThreeDFSC_Start.py')
        return str(cmd)


pyworkflow.em.Domain.registerPlugin(__name__)