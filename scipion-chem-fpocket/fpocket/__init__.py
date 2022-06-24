# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *           Lobna Ramadane Morchadi (lobna.ramadane@alumnos.upm.es)
# * Biocomputing Unit, CNB-CSIC
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

import pwem
from os.path import join, exists, abspath
from .constants import *

_version_ = '0.1'
_logo = "fpocket_logo.png"
_references = ['']


class Plugin(pwem.Plugin):
    _homeVar = FPOCKET_HOME
    _pathVars = [FPOCKET_HOME]
    _supportedVersions = [V3_0]
    _pluginHome = join(pwem.Config.EM_ROOT, FPOCKET + '-' + FPOCKET_DEFAULT_VERSION)

    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file.
        """
        cls._defineEmVar(FPOCKET_HOME, FPOCKET + '-' + FPOCKET_DEFAULT_VERSION)

    @classmethod
    def defineBinaries(cls, env):
        installationCmd = ''
        print('Installing with conda')
        installationCmd += 'git clone https://github.com/Discngine/fpocket.git {} && cd {} && make && '.format(cls._pluginHome, cls._pluginHome)

        # Creating validation file
        FPOCKET_INSTALLED = '%s_installed' % FPOCKET
        installationCmd += 'touch %s' % FPOCKET_INSTALLED  # Flag installation finished

        env.addPackage(FPOCKET,
                       version=FPOCKET_DEFAULT_VERSION,
                       tar='void.tgz',
                       commands=[(installationCmd, FPOCKET_INSTALLED)],
                       neededProgs=["conda"],
                       default=True)

    @classmethod
    def runFpocket(cls, protocol, program, args, cwd=None):
        """ Run Fpocket command from a given protocol. """
        protocol.runJob(join(cls._pluginHome, 'bin/{}'.format(program)), args, cwd=cwd)

    @classmethod
    def runMDpocket(cls, protocol, program, args, cwd=None):
        """ Run MDpocket command from a given protocol. """
        chg = 'mkdir {}/run_{} && cd {}/run_{} && {}'.\
            format(join(cls._pluginHome, 'bin'), protocol.getObjId(), join(cls._pluginHome, 'bin'), protocol.getObjId(),
                   program)
        print(args)
        protocol.runJob(chg, args, cwd=cwd)
        cmd = 'cd {}/run_{} && mv ./* {} && rm -r {}/run_{}'.format(join(cls._pluginHome, 'bin'), protocol.getObjId(),
                                                                    abspath(cwd),
                                                                    join(cls._pluginHome, 'bin'), protocol.getObjId())

        protocol.runJob(cmd, '', cwd=cwd)

    @classmethod
    def runSelIsovalue(cls, protocol, program, args, cwd=None):
        cmd = 'python {}/{}'.format(join(cls._pluginHome, 'scripts'), program)
        protocol.runJob(cmd, args, cwd=cwd)


    @classmethod
    def runMDpocket_2(cls, protocol, program, args, cwd=None ):
        protocol.runJob(program, args, cwd=cwd)

    @classmethod  #  Test that
    def getEnviron(cls):
        pass


    # ---------------------------------- Utils functions  -----------------------

