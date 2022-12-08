# corespray
A python package for sampling a distribution function for stars that have been ejected from a star cluster's core. If you use corespray in your research, please cite [Grondin et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.tmp.3150G/abstract) and link to https://github.com/webbjj/corespray/.

# Installation 
To install corespray from GitHub, clone the repository and install via setup tools:

git clone https://github.com/webbjj/corespray.git  
cd corespray  
python setup.py install  

Please note that if you don’t have permission to write files to the default install location (often /usr/local/lib), you will either need to run:

sudo python setup.py install 

or

python setup.py install --prefix='PATH' 

where ‘PATH’ is a directory that you do have permission to write in and is in your PYTHONPATH.

# Requirements  

corespray requires the following python packages:  

galpy (https://docs.galpy.org/en/v1.8.1/)  
matplotlib   
numpy  
