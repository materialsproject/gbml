from setuptools import setup, find_packages, Extension
import numpy
from codecs import open
from os import path

# References:
#   https://github.com/pypa/sampleproject/blob/master/setup.py

# Get the long description from the relevant file
here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'DESCRIPTION.txt'), encoding='utf-8') as file:
    long_description = file.read()

source_list=['src/band.c', 'src/dbinom.c', 'src/dens_haz.c', 'src/dens_int.c', 'src/density.c',
    'src/dens_odi.c', 'src/ev_atree.c', 'src/ev_interp.c', 'src/ev_kdtre.c', 'src/ev_main.c',
    'src/ev_sphere.c', 'src/ev_trian.c', 'src/family.c', 'src/fitted.c', 'src/frend.c',
    'src/gbml.c', 'src/lf_adap.c', 'src/lf_dercor.c', 'src/lf_fitfun.c', 'src/lf_nbhd.c',
    'src/lf_robust.c', 'src/lfstr.c', 'src/lf_vari.c', 'src/lf_wdiag.c', 'src/locfit.c',
    'src/math.c', 'src/m_chol.c', 'src/m_eigen.c', 'src/m_icirc.c', 'src/m_imont.c',
    'src/minmax.c', 'src/m_isimp.c', 'src/m_isphr.c', 'src/m_jacob.c', 'src/m_max.c',
    'src/m_qr.c', 'src/m_solve.c', 'src/m_svd.c', 'src/m_vector.c', 'src/pcomp.c',
    'src/predict.c', 'src/preplot.c', 'src/prob.c', 'src/procv.c', 'src/readGBML.c',
    'src/scb.c', 'src/scb_cons.c', 'src/scb_crit.c', 'src/scb_iface.c', 'src/simul.c',
    'src/smisc.c', 'src/spreplot.c', 'src/startlf.c', 'src/weight.c']

header_list=['src/cversion.h', 'src/design.h', 'src/imatlb.h', 'src/lfcons.h', 'src/lffuns.h',
    'src/lfstruc.h', 'src/lfwin.h', 'src/local.h', 'src/mutil.h', 'src/notR.h', 'src/predict.h',
    'src/readGBML.h', 'src/spreplot.h', 'src/tube.h']

core=Extension(
    'gbml.core',
    define_macros=[('MAJOR_VERSION', '1'), ('MINOR_VERSION', '0')],
    sources=source_list,
    depends=header_list,
    include_dirs=[numpy.get_include()])

setup(
    name='gbml',
    version='1.0.3',
    description='GBM-Locfit: A GBM framework using Locfit',
    long_description=long_description,
    author='Randy Notestine',
    author_email='RNotestine@UCSD.edu',
    license='GPL2',
    packages=find_packages(),
    url='https://github.com/materialsproject/gbml/',
    install_requires=['numpy', 'pymatgen'],
    include_package_data=True,
    ext_modules=[core],
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 2.7',
    ]
)

#   package_data={'gbml': ['gbml/data/*.data', 'gbml/data/*.json', 'src/*.h', 'src/README']},
