# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Lobna Ramadane Morchadi (lobna.ramadane@alumnos.upm.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************


"""
This protocol is used to perform a pocket search on a protein structure using the FPocket software

"""

import os, shutil

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pyworkflow.object import String
from pwem.protocols import EMProtocol
#from pwem.protocols import protocol_define_manual_pockets

from pwchem.objects import SetOfPockets, PredictPocketsOutput, ProteinPocket
from pwchem.utils import *
from fpocket import Plugin
from fpocket.constants import *

class MDpocketAnalyze(EMProtocol):
    """
    Executes the mdpocket software to look for protein pockets.
    """
    _label = 'Analyze pockets'
    _pocketTypes = ['Default Pockets', 'Druggable Pockets', 'Channels and small cavities', 'Water binding sites', 'Big external pockets']
    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSystem', params.PointerParam,
                       pointerClass='GromacsSystem', allowsNull=False,
                       label="Input atomic system: ",
                       help='Select the atom structure to search for pockets')

        form.addSection(label='Pocket analysis parameters')
        group = form.addGroup('Features')
        group.addParam('isoValue', params.FloatParam, default=1.0,
                       label='Selected Isovalue',
                       help='Selected Isovalue Threshold in Pocket Analysis for PDB output')
        group.addParam('maxIntraDistance', params.FloatParam, default='2.0',
                       label='Maximum distance between pocket points (A): ',
                       help='Maximum distance between two pocket atoms to considered them same pocket')
        group = form.addGroup('Pocket Type')
        group.addParam('pockType', params.EnumParam,
                   choices=self._pocketTypes, default=0,
                   label='Choose the type of pocket:',
                   help='Detect different type of pockets with a set of specific inner parameters'
                   )

    def _getMDpocketArgs(self):
        trajFile = os.path.abspath(self.inputSystem.get().getTrajectoryFile())
        args = ['--trajectory_file', trajFile]

        trajExt = os.path.splitext(trajFile)[1][1:]
        args += ['--trajectory_format', trajExt]

        pdbFile = os.path.abspath(self.inputSystem.get().getSystemFile())
        args += ['-f', pdbFile]

        selPock = self.getEnumText('pockType')



        if selPock == 'Default Pockets':
            pass

        elif selPock == 'Druggable Pockets':
            args += ['-S']

        elif selPock == 'Channels and small cavities':
            args += [' -m 2.8 -M 5.5 -i 3']

        elif selPock == 'Water binding sites':
            args += ['-m 3.5 -M 5.5 -i 3']

        elif selPock == 'Big external pockets':
            args += ['-m 3.5 -M 10.0 -i 3']

        args +=['-C']

        return args

    def _getselIsovalueArgs(self):
    # python extractISOPdb.py path / my_dx_file.dx outputname.pdb isovalue

        volFile = os.path.abspath(self._getExtraPath("mdpout_dens_grid.dx")) #To get the extra path between home path and protocol
        args = [volFile]

        outputName = 'mdpoutput-{}.pdb'.format(str(self.isoValue.get()))
        args += [outputName]

        isoValue = self.isoValue.get()
        args += [isoValue]
        return args

    def _clusterizedPocketsArgs(self):
        inpPdb = os.path.abspath(self._getExtraPath("mdpout_dens_iso_8.dx"))
        inpPdb += [inpPdb]
        return args

    #
    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('mdPocketStep')
        self._insertFunctionStep('createOutputStep')
        self._insertFunctionStep('selIsovalue')
        self._insertFunctionStep('defineOutputStep')

    def mdPocketStep(self):
        Plugin.runMDpocket(self, 'mdpocket', args=self._getMDpocketArgs(), cwd=self._getExtraPath())

    def selIsovalue(self):
        Plugin.runSelIsovalue(self, 'extractISOPdb.py', args=self._getselIsovalueArgs(), cwd=self._getExtraPath())

    def defineOutputStep(self):
        coords = self.getCoords()
        self.coordsClusters = clusterCoords(coords, self.maxIntraDistance.get())

        outPockets = SetOfPockets(filename=self._getPath('pockets.sqlite'))
        for i, clust in enumerate(self.coordsClusters):
            pocketFile = self.createPocketFile(clust, i)
            pocket = ProteinPocket(pocketFile, self.inputSystem.get().getSystemFile())
            pocket.calculateContacts()
            outPockets.append(pocket)

        if len(outPockets) >= 0: #Sometimes with the isovalue of 1 no pockets are detected, so still we want the output to be visualized
            outPockets.buildPDBhetatmFile()
            outPockets.densVolFile = String(os.path.abspath(self._getExtraPath('mdpout_dens_grid.dx')))
            outPockets.freqVolFile = String(os.path.abspath(self._getExtraPath('mdpout_freq_grid.dx')))
            self._defineOutputs(outputPockets=outPockets)

    def createOutputStep(self):
        pass


    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _warnings(self):
        """ Try to find warnings on define params. """
        warnings = []
        return warnings

    # --------------------------- UTILS functions -----------------------------------

    def getCoords(self):
        pdbFile = os.path.abspath(self._getExtraPath('mdpoutput-{}.pdb'.format(str(self.isoValue.get()))))
        coords = []
        for pdbLine in open(pdbFile):
            line = pdbLine.split()
            coord = line[5:8]
            coord = list(map(float,coord))
            coords.append(coord)
        return coords

    def createPocketFile(self, clust, i):
        outFile = self._getExtraPath('pocketFile_{}.pdb'.format(i+1))
        with open(outFile, 'w') as f:
            for j, coord in enumerate(clust):
                f.write(writePDBLine(['HETATM', str(j + 1), 'APOL', 'STP', 'C', '1', *coord, 1.0, 0.0, '', 'Ve']))
                #f.write(writePDBLine(['ATOM', str(j+1), '', 'C', 'PTH', '1', *coord, 0.0, 0.0, '', 'C']))

            f.write('\nEND')
        return outFile

    # def clusterCoords(coords, maxDist):
    #     clusters = []
    #     for coord in coords:
    #         newClusters = []
    #         newClust = [coord]
    #         for clust in clusters:
    #             merge = False
    #             for cCoord in clust:
    #                 dist = calculateDistance(coord, cCoord)
    #                 if dist < maxDist:
    #                     merge = True
    #                     break
    #
    #             if merge:
    #                 newClust += clust
    #             else:
    #                 newClusters.append(clust)
    #
    #         newClusters.append(newClust)
    #         clusters = newClusters.copy()
    #     return clusters