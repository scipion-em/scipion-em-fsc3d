# NYSBC - 3DFSC plugin

This plugin provide a wrapper around [3DFSC](https://github.com/nysbc/Anisotropy) program developed at Salk Institute and NYSBC.

![build status](http://scipion-test.cnb.csic.es:9980/badges/nysbc_devel.svg "Build status")

## Installation

You will need to use [2.0](https://github.com/I2PC/scipion/releases/tag/v2.0) version of Scipion to be able to run these protocols. To install the plugin, you have two options:

   a) Stable version
   ```
   scipion installp -p scipion-em-nysbc
   ```
   b) Developer's version
   * download repository 
   ```
    git clone https://github.com/scipion-em/scipion-em-nysbc.git
   ```
   * install 
   ```
    scipion installp -p path_to_scipion-em-nysbc --devel
   ```

3DFSC binaries will be installed automatically with the plugin, but you can also link an existing installation. **Important:** you need to have conda (miniconda3 or anaconda3) pre-installed to use this program. 
Default installation path assumed is `software/em/nysbc3DFSC-2.5`, if you want to change it, set *NYSBC_3DFSC_HOME* in `scipion.conf` file to the folder where the 3DFSC is installed. 

To check the installation, simply run the following Scipion test: `scipion test nysbc.tests.test_protocols_nysbc.TestNysbc3DFSC`

## Supported versions

2.5, 3.0

## Protocols

* 3d fsc

## References
1.  Yong Zi Tan, Philip R Baldwin, Joseph H Davis, James R Williamson, Clinton S Potter, Bridget Carragher & Dmitry Lyumkis. Addressing preferred specimen orientation in single-particle cryo-EM through tilting. Nature Methods volume 14, pages 793–796 (2017).
