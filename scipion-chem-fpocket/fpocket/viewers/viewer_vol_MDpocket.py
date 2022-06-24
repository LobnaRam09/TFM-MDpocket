# **************************************************************************
# *
# * Authors:     Lobna Ramadane Morchadi (lobna.rm@outlook.com)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

from ..protocols import MDpocketAnalyze
import pyworkflow.protocol.params as params
from pwchem.viewers import ViewerGeneralPockets
from pwchem.viewers import VmdViewPopen


class viewerVolMDPocket(ViewerGeneralPockets):
  _label = 'Viewer volumentric files of pockets MDPocket'
  _targets = [MDpocketAnalyze]
  _volFile = ['Frequency grid file', 'Density grid file']

  def __init__(self, **kwargs):
    super().__init__(**kwargs)

  def _defineParams(self, form):
    super()._defineParams(form)
    form.addSection(label='VMD visualization')
    form.addParam('displayVMD', params.LabelParam,
                  label='Display output pockets with VMD: ',
                  help='*VMD*: display output pockets and movies with VMD.')

    group = form.addGroup('Volumetric Files')
    group.addParam('displayVolFile', params.EnumParam,
                   choices=self._volFile, default=0,
                 label='Output volumetric files to visualize:',
                 help='*Frequency grid*: measure of how many times the pocket was open during MD trajectory\n'
                      '*Density grid*: Superposition of the alpha spheres of all snapshots along the MD trajectory'
                 )

  def _getVisualizeDict(self):
    visDict = super()._getVisualizeDict()
    visDict.update({'displayVMD': self._showVolFileVMD})
    return visDict

  def _validate(self):
    return []

  # =========================================================================

  def _showVolFileVMD(self, paramName=None):

    pdbFile = self.protocol.inputSystem.get().getSystemFile()
    #print('Dens file: ', self.protocol.outputPockets.densVolFile.get())
    densFile = self.protocol.outputPockets.densVolFile.get()
    freqFile = self.protocol.outputPockets.freqVolFile.get()


    TCL_MD_STR = '''
    mol representation Isosurface
    mol addrep 0
    mol new {%s} type {dx} first 0 last -1 step 1 waitfor 1 volsets {0 }
    animate style Loop
    menu files off
    menu files off
    menu files on
    display resetview
    mol addrep 1
    display resetview
    mol new {%s} type {pdb} first 0 last -1 step 1 waitfor 1 volsets {0 }
    animate style Loop
    menu files off
    menu graphics off
    menu graphics on
    mol modstyle 0 0 Isosurface 0.500000 0 2 1 1 1
    mol modstyle 0 1 NewCartoon 0.300000 10.000000 4.100000 0
    '''
    outTcl = self.protocol._getExtraPath('vmdSimulation.tcl')
    volFile = self.getEnumText('displayVolFile')

    if volFile == 'Frequency grid file':
        gridFile = freqFile
    elif volFile == 'Density grid file':
        gridFile = densFile

    with open(outTcl, 'w') as f:
      f.write(TCL_MD_STR % (gridFile, pdbFile))
    args = '-e {}'.format(outTcl)

    return [VmdViewPopen(args)]


