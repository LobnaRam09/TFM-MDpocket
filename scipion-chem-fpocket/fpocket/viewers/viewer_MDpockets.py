# **************************************************************************
# *
# * Authors:     Lobna Ramadane Morchadi (lobna.ramadane@alumnos.upm.es)
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
import os, glob

import matplotlib.pyplot as plt
import numpy as np
import pyworkflow.viewer as pwviewer


from ..protocols import MDpocketCharacterize
import pyworkflow.protocol.params as params
from pwchem.viewers import  PyMolView
from pwchem.utils import natural_sort
from pwchem.viewers import VmdViewPopen

from pwem.viewers.plotter import EmPlotter, plt
from matplotlib import ticker


from matplotlib.ticker import Formatter

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)


class viewerMDPocket(pwviewer.ProtocolViewer):
  _label = 'Viewer dynamic pockets of MDPocket'
  _targets = [MDpocketCharacterize]
  _descriptors = ['Pocket Volume', 'Solvent accessible surface area', 'Polar solvent accessible surface area', 'Apolar solvent accessible surface area',
                  'Accessible surface area (probe of 2.2A)', 'Polar solvent accessible surface area (probe of 2.2A)', 'Apolar solvent accessible surface area (probe of 2.2A)',
                  'Number of alpha-spheres', 'Mean alpha sphere radius', 'Mean alpha sphere solvent accessibility ','Pocket proportion of apolar alpha spheres', 'Mean local hydrophobic density',
                  'Hidrophobicity Score', 'Volume Score', 'Polarity Score', 'Charge Score', 'Proportion of polar atoms', 'Alpha sphere density', 'Max distance between mass center and all alpha spheres']
  def __init__(self, **kwargs):
      pwviewer.ProtocolViewer.__init__(self, **kwargs)

  def _defineParams(self, form):
    form.addSection(label='Pymol visualization')
    group = form.addGroup('Pockets to view')
    group.addParam('nPocket', params.EnumParam, default=0, #enumParam para ver la lista de carpetas
                   choices= self._getDynPockets(),
                   label='Choose the pocket to visualize:',
                   help='Selected pockets and protein atom interactions to visualize along the MD trajectory')


    group = form.addGroup('Open MD simulation')

    group.addParam('displayPymol', params.LabelParam,
                  label='Display dynamic pocket with Pymol: ',
                  help='*Pymol*: display dynamic Pockets along the MD simulation.')

    group.addParam('displayMdVMD', params.LabelParam,
                   label='Display receptor atoms of pocket with VMD: ',
                   help='Display receptor atoms defining the binding pocket interaction in MD trajectory with VMD.'

                   )


    group = form.addGroup('Pockets Description')
    group.addParam('displayGraphic', params.LabelParam,
                  label='Display a graph of selected pocket: ',
                  help='Display a graphical representation of descriptors of selected pockets along MD trajectory')
    group.addParam('displayDesc', params.EnumParam,
                   choices=self._descriptors, default=0,
                   label='Choose the descriptor to visualize:',
                   help='Graphical representation of descriptor variation along the MD trajectory of previous selected pocket'
                   )
  def _getVisualizeDict(self):
      return {
        'displayPymol': self._showMdPymol,
        'displayMdVMD': self._showMdVMD,
        'displayDesc': self._displayGraphic
      }

  def _validate(self):
    return []

  # =========================================================================

  def _showMdPymol(self, paramName=None):
    PYMOL_MD_POCK = ''' load {}
    load_traj {}
    set movie_fps, 15
    load {}
    load {}
    
    '''
    pdbFile = self.protocol.inputSystem.get().getSystemFile()
    trjFile = self.protocol.inputSystem.get().getTrajectoryFile()
    dir = os.path.abspath(self.protocol._getExtraPath('pocketFolder_{}'.format(str(self.nPocket.get()+1))))
    dynPocket ='{}/mdpout_mdpocket_{}.pdb'.format(dir, str(self.nPocket.get()+1))
    dynAtoms = '{}/mdpout_mdpocket_atoms_{}.pdb'.format(dir, str(self.nPocket.get()+1))
    outPml = self.protocol._getExtraPath('pymolSimulation.pml')
    with open(outPml, 'w') as f:
      f.write(PYMOL_MD_POCK .format(os.path.abspath(pdbFile),
                                    os.path.abspath(trjFile),
                                    os.path.abspath(dynPocket),
                                    os.path.abspath(dynAtoms)))



    return [PyMolView(os.path.abspath(outPml), cwd=self.protocol._getPath())]


  def _getDynPockets(self):
      n_pockets = []
      for pockDir in natural_sort(glob.glob(self.protocol._getExtraPath('pocketFolder_*'))):
          n_pocket = os.path.basename(pockDir)
          n_pocket = n_pocket.split('_')[1]
          n_pockets.append(n_pocket)

      return n_pockets



  def _showMdVMD(self, paramName=None):
      dir = os.path.abspath(self.protocol._getExtraPath('pocketFolder_{}'.format(str(self.nPocket.get() + 1))))
      pdbFile = self.protocol.inputSystem.get().getSystemFile()
      dynAtoms = '{}/mdpout_mdpocket_atoms_{}.pdb'.format(dir, str(self.nPocket.get() + 1))
      trjFile = self.protocol.inputSystem.get().getTrajectoryFile()

      TCL_MD_STR = '''
        mol addrep 0
        mol new {%s} type {pdb} first 0 last -1 step 1 waitfor 1
        mol addfile {%s} type {xtc} first 0 last -1 step 1 waitfor 1 0
        
        mol color Name
        mol representation NewCartoon 0.300000 10.000000 4.100000 0
        mol selection protein
        mol material Opaque
        mol modrep 0 0
        
        mol addrep 0
        mol color Name
        mol representation Points 1.000000
        mol selection hetero within 3 of protein
        mol material Opaque
        mol modrep 1 0
        
        
        display resetview
        mol new {%s} type {pdb} first 0 last -1 step 1 waitfor 1 0
        mol modstyle 0 1 Licorice 0.300000 12.000000 12.000000
        
        '''

      outTcl = self.protocol._getExtraPath('vmdSimulation.tcl')
      with open(outTcl, 'w') as f:
          f.write(TCL_MD_STR % (pdbFile, trjFile, dynAtoms))
      args = '-e {}'.format(outTcl)

      return [VmdViewPopen(args)]



  def _displayGraphic(self,  paramName=None):

      dir = os.path.abspath(self.protocol._getExtraPath('pocketFolder_{}'.format(str(self.nPocket.get()+1))))
      descrFile = '{}/mdpout_descriptors_{}.txt'.format(dir, str(self.nPocket.get()+1))

      with open(descrFile, 'r') as fileTxt:
          desc_dic = {} #dictionary descriptor:value
          header = fileTxt.readline().split() #create a list with the header of the file
          for h in header:
              desc_dic[h] = [] # inicializate the dictionary so we can append the values
          for line in fileTxt:
              for i, descriptor in enumerate(line.split()): #take the i (position) and descriptor value in each line of the file
                  desc_dic[header[i]].append(descriptor) #add to the dictionary the value of each descriptor




      desc_dic_list = list(desc_dic.values())
      snaps = desc_dic_list[0] #Snapshots are the first element of the list
      d = self.displayDesc.get()+1 #To start by the second element
      desctr = list(map(float, desc_dic_list[d]))



      # Plot of the descriptors vs snapsohts:
      self.plotter = EmPlotter(x = 1, y = 1, windowTitle='Pocket Descriptors')
      a = self.plotter.createSubPlot("Pocket {} ".format(str(self.nPocket.get()+1)), "Snapshots", "Descriptors")

      a.plot(snaps, desctr)

      a.legend(labels=['{}'.format(str(self.getEnumText('displayDesc')))]) #Get the name of the list element  displayed

      # Formatting the X axis and Y axis for correct values distribution in the axis
      xmajor_formatter = FormatStrFormatter('%1.1f') # 1 space reserved for decimal value
      ymajor_formatter = FormatStrFormatter('%1.1f')

      a.xaxis.set_major_formatter(xmajor_formatter)

      a.yaxis.set_major_formatter(ymajor_formatter)

      a.yaxis.set_major_locator(ticker.AutoLocator())
      #a.yaxis.set_minor_locator(ticker.AutoMinorLocator())


      if len(snaps) < 500:
          a.xaxis.set_major_locator(ticker.IndexLocator(offset=0, base=5))


      else:
          locator = ticker.AutoLocator()
          locator.set_params(integer=True)
          minor_locator = ticker.AutoMinorLocator()

          a.xaxis.set_major_locator(locator)
          a.xaxis.set_minor_locator(minor_locator)




      return [self.plotter]


