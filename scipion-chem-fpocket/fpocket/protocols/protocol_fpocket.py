# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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

from pwchem.objects import SetOfPockets, PredictPocketsOutput, ProteinPocket
from pwchem.utils import clean_PDB

from fpocket import Plugin
from fpocket.constants import *

class FpocketFindPockets(EMProtocol):
    """
    Executes the fpocket software to look for protein pockets.
    """
    _label = 'Find pockets'
    _possibleOutputs = PredictPocketsOutput

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputAtomStruct', params.PointerParam,
                       pointerClass='AtomStruct', allowsNull=False,
                       label="Input atom structure",
                       help='Select the atom structure to search for pockets')

        form.addSection(label='Pocket detection parameters')
        group = form.addGroup('Alpha spheres')
        group.addParam('minAlpha', params.FloatParam, default=3.4,
                       label='Min alpha sphere radius',
                       help='Minimum radius of an alpha-sphere (A)')
        group.addParam('maxAlpha', params.FloatParam, default=6.2,
                      label='Max alpha sphere radius',
                      help='Maximum radius of an alpha-sphere (A)')
        group.addParam('minNSpheres', params.IntParam, default=15,
                       label='Min a-spheres per pocket',
                       help='Minimum number of a-sphere per pocket')
        group.addParam('ratioApSpheres', params.FloatParam, default=0.0,
                       label='Min ratio of a-spheres', expertLevel=params.LEVEL_ADVANCED,
                       help='Minimum proportion of apolar sphere in a pocket (remove otherwise)')
        group.addParam('minApNeigh', params.IntParam, default=3,
                       label='Min apolar neigh per sphere', expertLevel=params.LEVEL_ADVANCED,
                       help='Minimum number of apolar neighbor for an a-sphere to be considered as apolar.')
        group = form.addGroup('Clustering')
        group.addParam('clustType', params.EnumParam,
                       choices=CLUST_TYPES, label="Clustering linkage type", default=0,
                       help='Specify the clustering method wanted for grouping voronoi vertices together')
        group.addParam('clustDistType', params.EnumParam,
                       choices=DIST_TYPES, label="Clustering distance type", default=0,
                       help='Specify the distance measure for clustering')
        group.addParam('clustDist', params.FloatParam, default=2.4,
                       label='Clustering distance threshold',
                       help='Distance threshold for clustering algorithm')
        form.addParam('mcIterVol', params.IntParam, default=300, expertLevel=params.LEVEL_ADVANCED,
                       label='Monte-Carlo iterations for volume',
                       help='Number of Monte-Carlo iteration for the calculation of each pocket volume.')


    def _getFpocketArgs(self):
        args = ['-f', os.path.abspath(self.inpFile)]

        #Alpha spheres
        args += ['-m', self.minAlpha.get(), '-M', self.maxAlpha.get(), '-i', self.minNSpheres.get(),
                 '-p', self.ratioApSpheres.get(), '-A', self.minApNeigh.get()]
        #Clustering
        args += ['-C', CLUST_TYPES_CODES[self.clustType.get()],
                 '-e', DIST_TYPES_CODES[self.clustDistType.get()],
                 '-D', self.clustDist.get()]
        #Volume
        args += ['-v', self.mcIterVol.get()]

        return args

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('fPocketStep')
        self._insertFunctionStep('createOutputStep')

    def convertInputStep(self):
        #Simply copying the input struct file into current extra file (fpocket will create output there automatically)
        inpStruct = self.inputAtomStruct.get()
        inpFile = self.getPdbInputStruct()
        inpName = self.getPdbInputStructName()
        self.inpBase, self.ext = os.path.splitext(inpName)
        if self.ext == '.ent':
            self.inpFile = self._getExtraPath(self.inpBase+'.pdb')
            shutil.copy(inpFile, self.inpFile)
        elif self.ext == '.cif':
            self.inpFile = self._getExtraPath(self.inpBase + '.pdb')
            clean_PDB(inpFile, self.inpFile)
        elif str(type(inpStruct).__name__) == 'SchrodingerAtomStruct':
            self.inpFile = self._getExtraPath(os.path.basename(inpFile).replace(inpStruct.getExtension(), '.pdb'))
            inpStruct.convert2PDB(outPDB=self.inpFile)
        else:
            self.inpFile = self._getExtraPath(inpName)
            shutil.copy(inpFile, self.inpFile)

    def fPocketStep(self):
        Plugin.runFpocket(self, 'fpocket', args=self._getFpocketArgs(), cwd=self._getExtraPath())

    def createOutputStep(self):
        pocketsDir = os.path.abspath(self._getExtraPath('{}/pockets'.format(self.inpBase+'_out')))
        pocketFiles = os.listdir(pocketsDir)

        inpStruct = self.inputAtomStruct.get()
        outPockets = SetOfPockets(filename=self._getExtraPath('pockets.sqlite'))
        for pFile in pocketFiles:
            if '.pdb' in pFile:
                pFileName = os.path.join(pocketsDir, pFile)
                pqrFile = pFileName.replace('atm.pdb', 'vert.pqr')
                pock = ProteinPocket(pqrFile, self.inpFile, pFileName, pClass='FPocket')
                if str(type(inpStruct).__name__) == 'SchrodingerAtomStruct':
                  pock._maeFile = String(os.path.abspath(inpStruct.getFileName()))
                outPockets.append(pock)

        outHETMFile = outPockets.buildPDBhetatmFile()
        self._defineOutputs(**{self._possibleOutputs.outputPockets.name: outPockets})


    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _warnings(self):
        """ Try to find warnings on define params. """
        import re
        warnings = []
        inpStruct = self.inputAtomStruct.get()
        inpFile = os.path.abspath(inpStruct.getFileName())
        if str(type(inpStruct).__name__) == 'SchrodingerAtomStruct':
          inpFile = inpStruct.convert2PDB()
        with open(inpFile) as f:
          fileStr = f.read()
        if re.search('\nHETATM', fileStr):
          warnings.append('The structure you are inputing has some *heteroatoms* (ligands).\n'
                          'This will affect the results as its volume is also taken as target.')

        return warnings

    # --------------------------- UTILS functions -----------------------------------
    def getPdbInputStruct(self):
        return self.inputAtomStruct.get().getFileName()

    def getPdbInputStructName(self):
        return self.getPdbInputStruct().split('/')[-1]

    def getPDBName(self):
        return self.getPdbInputStructName().split('.')[0]