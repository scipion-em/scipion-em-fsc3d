# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

from pyworkflow.protocol.params import LabelParam, EnumParam
from pyworkflow.viewer import DESKTOP_TKINTER
from pwem.viewers import ChimeraView, ObjectView, EmProtocolViewer

from .protocols import Prot3DFSC
from .constants import (VOLUME_SLICES, VOLUME_CHIMERA,
                        VOL_ORIG, VOL_TH, VOL_THBIN)


class ThreedFscViewer(EmProtocolViewer):
    """ Visualization of 3D FSC results. """
           
    _environments = [DESKTOP_TKINTER]
    _targets = [Prot3DFSC]
    _label = 'viewer'

    def __init__(self, **kwargs):
        EmProtocolViewer.__init__(self, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        group = form.addGroup('Volumes')
        group.addParam('displayVol', EnumParam, choices=['slices', 'chimera'],
                       default=VOLUME_SLICES, display=EnumParam.DISPLAY_HLIST,
                       label='Display volume with',
                       help='*slices*: display volumes as 2D slices along z axis.\n'
                            '*chimera*: display volumes as surface with Chimera.')
        group.addParam('doShowOutVol', EnumParam, default=VOL_ORIG,
                       choices=['original', 'thresholded',
                                'thresholded and binarized', 'all'],
                       display=EnumParam.DISPLAY_COMBO,
                       label='3D FSC volume to display')

        form.addParam('doShowHistogram', LabelParam,
                      label="Show histogram and directional FSC plot")
        form.addParam('doShowPlotFT', LabelParam,
                      label="Show Fourier Transform Power plot")
        form.addParam('doShowPlot3DFSC', LabelParam,
                      label="Show 3D FSC plots")
        form.addParam('doShowChimera', LabelParam,
                      label="Show Chimera animation", default=True,
                      help="Display 3D FSC and coloring original map by "
                           "angular resolution.")
        
    def _getVisualizeDict(self):
        self.protocol._initialize()  # Load filename templates
        return {'doShowOutVol': self._showVolumes,
                'doShowHistogram': self._showHistogram,
                'doShowPlotFT': self._showPlotFT,
                'doShowPlot3DFSC': self._showPlot3DFSC,
                'doShowChimera': self._showChimera
                }

# =============================================================================
# ShowVolumes
# =============================================================================
    def _showVolumes(self, paramName=None):
        if self.displayVol == VOLUME_CHIMERA:
            return self._showVolumesChimera()
        elif self.displayVol == VOLUME_SLICES:
            return self._createVolumesSqlite()

    def _showVolumesChimera(self):
        """ Create a chimera script to visualize selected volumes. """
        volumes = self._getVolumeNames()
        cmdFile = self.protocol._getExtraPath('chimera_volumes.cxc')
        with open(cmdFile, 'w+') as f:
            for vol in volumes:
                # We assume that the chimera script will be generated
                # at the same folder as 3DFSC volumes
                if os.path.exists(vol):
                    localVol = os.path.relpath(vol,
                                               self.protocol._getExtraPath())
                    f.write("open %s\n" % localVol)
            f.write('tile\n')
        view = ChimeraView(cmdFile)
        return [view]

    def _createVolumesSqlite(self):
        """ Write an sqlite with all volumes selected for visualization. """
        path = self.protocol._getExtraPath('3DFSC_viewer_volumes.sqlite')
        samplingRate = self.protocol.inputVolume.get().getSamplingRate()

        vols = self._getVolumeNames()
        files = []
        for vol in vols:
            if os.path.exists(vol):
                files.append(vol)
        self.createVolumesSqlite(files, path, samplingRate)
        return [ObjectView(self._project, self.protocol.strId(), path)]

# =============================================================================
    def _showPlot(self, fn):
        img = mpimg.imread(self.protocol._getFileName(fn))
        imgplot = plt.imshow(img)
        plt.axis('off')
        plt.show()
        return [imgplot]

    def _showHistogram(self, param=None):
        self._showPlot('out_histogram')

    def _showPlotFT(self, param=None):
        self._showPlot('out_plotFT')

    def _showPlot3DFSC(self, param=None):
        self._showPlot('out_plot3DFSC')

    def _showChimera(self, param=None):
        return [self.errorMessage('ChimeraX is not supported for this animation yet.',
                                  title="Visualization error")]
        # cmdFile = self.protocol._getFileName('out_cmdChimera')
        # view = ChimeraView(cmdFile)
        # return [view]

    def _getVolumeNames(self):
        volsFn = ['out_vol3DFSC', 'out_vol3DFSC-th', 'out_vol3DFSC-thbin']

        if self.doShowOutVol.get() == VOL_ORIG:
            vols = [self.protocol._getFileName(volsFn[0])]
        elif self.doShowOutVol.get() == VOL_TH:
            vols = [self.protocol._getFileName(volsFn[1])]
        elif self.doShowOutVol.get() == VOL_THBIN:
            vols = [self.protocol._getFileName(volsFn[2])]
        else:
            vols = [self.protocol._getFileName(f) for f in volsFn]

        return vols
