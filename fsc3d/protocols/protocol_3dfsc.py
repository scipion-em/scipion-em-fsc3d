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
from enum import Enum

import pyworkflow.protocol.params as params
from pyworkflow.constants import PROD
from pyworkflow.utils import cleanPath
from pwem.protocols import ProtAnalysis3D
from pwem.emlib.image import ImageHandler
from pwem.objects import Volume

from .. import Plugin


class outputs(Enum):
    outputVolume = Volume


class Prot3DFSC(ProtAnalysis3D):
    """
    Calculates directional resolution anisotropy in three-dimensional
    cryo-EM reconstructions using 3D Fourier Shell Correlation analysis.
    The protocol evaluates how resolution varies across different spatial
    orientations, allowing users to identify anisotropic signal quality
    and directional reconstruction biases. More info:
    https://github.com/nysbc/Anisotropy

    AI Generated:

    3D FSC Protocol (Prot3DFSC) — User Manual
        Overview

        The 3D FSC protocol evaluates directional resolution in cryo-EM
        reconstructions through the calculation of three-dimensional Fourier
        Shell Correlation volumes. Its primary purpose is to determine whether
        the resolution of a reconstructed map is isotropic or varies depending
        on spatial orientation. This analysis is biologically important because
        anisotropic resolution often reflects preferred particle orientations,
        incomplete angular sampling, specimen flexibility, or limitations in
        data acquisition and processing.

        In standard global FSC analysis, a single resolution estimate is
        reported for the entire reconstruction. However, many biological
        datasets contain directional differences in signal quality that are not
        captured by a global value alone. The 3D FSC protocol provides a more
        realistic interpretation of map quality by identifying regions or
        directions where the reconstruction may be weaker or less reliable.

        Inputs and General Workflow

        The protocol requires a full reconstructed volume together with two
        corresponding half maps. These half maps are essential because the
        directional FSC calculation relies on comparing independently refined
        reconstructions. The full map is used as the reference reconstruction
        for visualization and interpretation of the directional resolution
        results.

        Users may either provide the half maps explicitly or rely on half maps
        already associated with the input reconstruction. This flexibility is
        useful when working with reconstructions imported from external cryo-EM
        software packages or when the metadata already preserves half-map
        relationships internally.

        Optionally, a mask can be applied during the analysis. The mask defines
        which structural regions contribute to the directional FSC calculation.
        This is particularly important for large macromolecular assemblies,
        membrane proteins, or flexible complexes where solvent regions or mobile
        domains could otherwise dominate the anisotropy measurements.

        Biological Meaning of Directional Resolution

        Directional resolution analysis is especially valuable in cryo-EM
        studies where preferred orientation is suspected. When particles adopt
        limited angular distributions on the grid, some viewing directions are
        oversampled while others remain poorly represented. The resulting maps
        may therefore contain high resolution information along certain
        directions but weaker signal in others.

        For biological interpretation, anisotropy can affect the visibility of
        secondary structure features, ligand densities, membrane boundaries, or
        flexible domains. Understanding these directional limitations helps
        users assess the reliability of structural conclusions and avoid
        overinterpretation of weakly resolved regions.

        The protocol also provides sphericity measurements, which summarize how
        isotropic the directional FSC distribution is. Values closer to perfect
        spherical symmetry indicate more uniform resolution in all directions,
        whereas lower sphericity values suggest significant anisotropy.

        Cone Sampling and Angular Resolution

        Directional FSC calculations are performed by sampling Fourier space
        within angular cones. The cone angle controls the balance between
        directional sensitivity and statistical robustness. Smaller cone angles
        provide finer directional discrimination but may increase noise and
        instability. Larger cone angles produce smoother and more stable
        estimates at the cost of reduced directional specificity.

        In biological workflows, moderate cone angles are generally appropriate
        for routine assessment of anisotropy. Extremely small values are mainly
        useful for specialized investigations of highly directional datasets.

        FSC Thresholds and Sphericity Analysis

        The protocol allows users to define the FSC cutoff criterion used for
        directional resolution estimation. The commonly used 0.143 criterion is
        provided as the standard default because it is widely accepted in cryo-
        EM resolution validation workflows.

        Sphericity thresholds define which portions of the directional FSC
        volume contribute to anisotropy calculations. Evaluating multiple
        thresholds can reveal whether anisotropy changes across spatial
        frequencies. This may provide insight into overfitting, orientation
        assignment problems, or frequency-dependent reconstruction artifacts.

        High-pass filtering can additionally stabilize thresholding behavior by
        suppressing low-frequency fluctuations that may distort anisotropy
        measurements. This is particularly useful for noisy datasets or maps
        containing strong low-resolution background variations.

        Masking Strategies

        Mask selection strongly influences the biological relevance of the
        anisotropy analysis. Broad masks including excessive solvent regions
        may artificially increase isotropy by diluting structural signal,
        whereas overly restrictive masks can exaggerate anisotropy by focusing
        only on compact regions.

        For flexible assemblies, masks centered on the stable structural core
        often produce the most interpretable directional FSC measurements.
        Membrane proteins, elongated assemblies, and filamentous complexes may
        particularly benefit from carefully designed masks.

        Outputs and Interpretation

        After completion, the protocol produces a three-dimensional FSC volume
        together with visualization plots and directional resolution analyses.
        These outputs help users identify preferred directions of signal loss,
        assess reconstruction quality, and communicate anisotropy properties in
        publications or validation reports.

        The generated FSC volume can be visualized in molecular graphics
        software to inspect anisotropy geometrically. Thresholded outputs are
        particularly useful for observing whether the reconstruction retains
        spherical symmetry or displays elongated or flattened resolution
        distributions.

        Practical Recommendations

        In routine cryo-EM workflows, users should first examine the global FSC
        together with the directional FSC plots. Strong anisotropy may indicate
        preferred particle orientation, inaccurate alignment, insufficient tilt
        coverage, or specimen flexibility. Applying an appropriate mask often
        improves the interpretability of the results.

        When analyzing membrane proteins or elongated assemblies, anisotropy is
        common and should not automatically be interpreted as processing
        failure. Instead, the directional FSC analysis should be considered as
        complementary information that contextualizes the biological reliability
        of the reconstruction.

        For publication-quality validation, users are encouraged to inspect the
        directional FSC volume visually and compare the measured anisotropy
        against known characteristics of the specimen and acquisition geometry.

        Final Perspective

        The 3D FSC protocol provides an advanced framework for evaluating
        directional resolution in cryo-EM reconstructions. Rather than relying
        solely on a single global resolution estimate, it enables a more
        realistic understanding of how structural information is distributed
        throughout Fourier space. Careful interpretation of anisotropy,
        combined with appropriate masking and threshold selection, is essential
        for producing biologically meaningful conclusions from cryo-EM maps.
    """
    _label = 'estimate resolution'
    _devStatus = PROD
    _possibleOutputs = outputs

    INPUT_HELP = """ Required input volumes for 3D FSC:
        1. First half map of 3D reconstruction. Can be masked or unmasked.
        2. Second half map of 3D reconstruction. Can be masked or unmasked.
        3. Full map of 3D reconstruction. Can be masked or unmasked, sharpened or unsharpened.
    """
    
    def __init__(self, **kwargs):
        ProtAnalysis3D.__init__(self, **kwargs)

    def _initialize(self):
        """ This function is mean to be called after the
        working dir for the protocol have been set. (maybe after recovery from mapper)
        """
        self._createFilenameTemplates()

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
                  'input_volFn': self._getTmpPath('volume_full.mrc'),
                  'input_half1Fn': self._getTmpPath('volume_half1.mrc'),
                  'input_half2Fn': self._getTmpPath('volume_half2.mrc'),
                  'input_maskFn': self._getTmpPath('mask.mrc'),
                  'out_histogram': self._getExtraPath('Results_vol/histogram.png'),
                  'out_plot3DFSC': self._getExtraPath('Results_vol/Plotsvol.jpg'),
                  'out_plotFT': self._getExtraPath('Results_vol/FTPlotvol.jpg'),
                  'out_vol3DFSC': self._getExtraPath('Results_vol/vol.mrc'),
                  'out_vol3DFSC-th': self._getExtraPath('Results_vol/vol_Thresholded.mrc'),
                  'out_vol3DFSC-thbin': self._getExtraPath('Results_vol/vol_ThresholdedBinarized.mrc'),
                  'out_cmdChimera': self._getExtraPath('Results_vol/Chimera/3DFSCPlot_Chimera.cmd'),
                  'out_globalFSC': self._getExtraPath('Results_vol/ResEMvolOutglobalFSC.csv')
                  }

        self._updateFilenamesDict(myDict)

    # --------------------------- DEFINE param functions ----------------------

    def _defineParams(self, form):
        form.addHidden(params.USE_GPU, params.BooleanParam,
                       default=True,
                       label="Use GPU?")
        form.addHidden(params.GPU_LIST, params.StringParam,
                       default='0',
                       label="Choose GPU ID",
                       help="Each GPU has a unique ID. If you have only "
                            "one GPU, set ID to 0. 3DFSC can use only one GPU.")

        form.addSection(label='Input')
        form.addParam('inputVolume', params.PointerParam,
                      pointerClass='Volume',
                      label="Input volume", important=True,
                      help=self.INPUT_HELP)
        form.addParam('provideHalfMaps', params.BooleanParam,
                      default=False,
                      label="Provide half-maps separately?")
        form.addParam('volumeHalf1', params.PointerParam,
                      label="Volume half 1", important=True,
                      condition="provideHalfMaps",
                      pointerClass='Volume',
                      help=self.INPUT_HELP)
        form.addParam('volumeHalf2', params.PointerParam,
                      pointerClass='Volume',
                      condition="provideHalfMaps",
                      label="Volume half 2", important=True,
                      help=self.INPUT_HELP)

        form.addParam('applyMask', params.BooleanParam, default=False,
                      label="Mask input volume?",
                      help='If given, it would be used to mask the half maps '
                           'during 3DFSC generation and analysis.')
        form.addParam('maskVolume', params.PointerParam, label="Mask volume",
                      pointerClass='VolumeMask', condition="applyMask",
                      help='Select a volume to apply as a mask.')

        form.addSection(label='Extra params')
        form.addParam('dTheta', params.FloatParam, default=20,
                      label='Angle of cone (deg)',
                      help='Angle of cone to be used for 3D FSC sampling in '
                           'degrees. Default is 20 degrees.')
        form.addParam('fscCutoff', params.FloatParam, default=0.143,
                      label='FSC cutoff',
                      help='FSC cutoff criterion. 0.143 is default.')
        form.addParam('thrSph', params.FloatParam, default=0.5,
                      label='Sphericity threshold',
                      help='Threshold value for 3DFSC volume for calculating '
                           'sphericity. 0.5 is default.')
        form.addParam('hpFilter', params.FloatParam, default=200,
                      label='High-pass filter (A)',
                      help='High-pass filter for thresholding in Angstrom. '
                           'Prevents small dips in directional FSCs at low '
                           'spatial frequency due to noise from messing up '
                           'the thresholding step. Decrease if you see a '
                           'huge wedge missing from your thresholded 3DFSC '
                           'volume. 200 Angstroms is default.')
        form.addParam('numThr', params.IntParam, default=0,
                      label='Number of threshold for sphericity',
                      help='Calculate sphericities at different threshold '
                           'cutoffs to determine sphericity deviation across '
                           'spatial frequencies. This can be useful to '
                           'evaluate possible effects of overfitting or '
                           'improperly assigned orientations. 0 is default.')

    # --------------------------- INSERT steps functions ----------------------
    
    def _insertAllSteps(self):
        # Insert processing steps
        self._initialize()
        self._insertFunctionStep('convertInputStep', needsGPU=False)
        self._insertFunctionStep('run3DFSCStep', needsGPU=self.usesGpu())
        self._insertFunctionStep('createOutputStep', needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    
    def convertInputStep(self):
        """ Convert input volumes to .mrc as expected by 3DFSC."""
        ih = ImageHandler()
        if self.provideHalfMaps:
            fnHalf1 = self.volumeHalf1.get().getLocation()
            fnHalf2 = self.volumeHalf2.get().getLocation()
        else:
            fnHalf1, fnHalf2 = self.inputVolume.get().getHalfMaps().split(',')

        ih.convert(fnHalf1,
                   self._getFileName('input_half1Fn'))
        ih.convert(fnHalf2,
                   self._getFileName('input_half2Fn'))
        ih.convert(self.inputVolume.get().getLocation(),
                   self._getFileName('input_volFn'))
        if self.maskVolume.hasValue():
            ih.convert(self.maskVolume.get().getLocation(),
                       self._getFileName('input_maskFn'))

    def run3DFSCStep(self):
        args = self._getArgs()
        params = ' '.join(['%s=%s' % (k, str(v)) for k, v in args.items()])

        if self.useGpu:
            params += ' --gpu --gpu_id=%s' % self.gpuList.get()

        Plugin.runProgram(self, params, cwd=self._getExtraPath())
        if not os.path.exists(self._getFileName('out_vol3DFSC')):
            raise RuntimeError('3D FSC run failed!')

    def createOutputStep(self):
        if os.path.exists(self._getFileName('out_vol3DFSC')):
            inputVol = self.inputVolume.get()
            vol = Volume()
            vol.setObjLabel('3D FSC')
            vol.setFileName(self._getFileName('out_vol3DFSC'))
            vol.setSamplingRate(inputVol.getSamplingRate())

            # remove useless output
            cleanPath(self._getExtraPath('Results_vol/ResEMvolOut.mrc'))

            self._defineOutputs(**{outputs.outputVolume.name: vol})
            self._defineSourceRelation(self.inputVolume, vol)

    # --------------------------- INFO functions ------------------------------
    
    def _summary(self):
        summary = []
        if self.getOutputsSize() > 0:
            logFn = self.getLogPaths()[0]
            sph = self.findSphericity(logFn)
            summary.append(f'Sphericity: {sph:0.3f}')
        else:
            summary.append("Output is not ready yet.")

        return summary
    
    def _validate(self):
        errors = []

        if not self.provideHalfMaps and not self.inputVolume.get().hasHalfMaps():
            errors.append("Input volume has no associated half-maps.")
                
        return errors
    
    # --------------------------- UTILS functions -----------------------------
 
    def _getArgs(self):
        """ Prepare the args dictionary."""

        args = {'--halfmap1': os.path.relpath(self._getFileName('input_half1Fn'),
                                              self._getExtraPath()),
                '--halfmap2': os.path.relpath(self._getFileName('input_half2Fn'),
                                              self._getExtraPath()),
                '--fullmap': os.path.relpath(self._getFileName('input_volFn'),
                                             self._getExtraPath()),
                '--apix': self.inputVolume.get().getSamplingRate(),
                '--ThreeDFSC': 'vol',
                '--dthetaInDegrees': self.dTheta.get(),
                '--FSCCutoff': self.fscCutoff.get(),
                '--ThresholdForSphericity': self.thrSph.get(),
                '--HighPassFilter': self.hpFilter.get(),
                '--numThresholdsForSphericityCalcs': self.numThr.get(),
                '--histogram': 'histogram'
                }
        if self.applyMask and self.maskVolume:
            args['--mask'] = os.path.relpath(self._getFileName('input_maskFn'),
                                             self._getExtraPath())
        return args

    def findSphericity(self, fn):
        with open(fn, 'r') as f:
            sph = 0.
            for line in f:
                if 'Sphericity is ' in line:
                    sph = float(line.split()[2])

        return sph
