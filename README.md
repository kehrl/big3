# big3
Scripts and utilities for running, processing, and plotting satellite observations and finite element model results for Helheim and Kangerlussuaq Glaciers, SE Greenland.

Code developed as part of my dissertation at the University of Washington.

## Dependencies
- Python (numpy, scipy, matplotlib, shapely, gdal, paraview, netcdf, and many others)
- [Elmer/Ice parallel version](https://github.com/ElmerCSC/elmerfem) (with fortran, open-mpi, mumps, hypre) for big3/modeling
- [Gmsh](http://gmsh.info/) for big3/modeling 


## To run code
- `big3/lib` includes all python libraries and will need to be bundled or added to python path
- elmer functions and modules in `big3/modeling/elmerlib` will need to be compiled and added to `DYLD_LIBRARY_PATH` to run simulations in Elmer 
- paths are presently hardcoded, so you will need to have a similar data and modeling directory structure to automatically load data and save model outputs (environment variables `MODEL_HOME`, `DATA_HOME`, and `CODE_HOME` will need to be set)

## big3/lib
- `bedlib.py`,`climlib.py`,`icefrontlib.py`,`vellib.py`, and `zslib.py` - pull, process, and output bed elevations, climate data (RACMO, OSTIA SIF/SST), ice-front positions, velocity measurements, and surface-elevation measurements
- `coordlib.py` - convert coordinates and get geoid heights
- `crossoverlib.py` - find crossovers and their differences in an array of three-dimensional points (x,y,z)
- `datelib.py` - convert between fractional year, date, and day of year
- `distlib.py` - find distance along a transect or between points
- `elmerreadlib.py` - read, process, and output Elmer model results from saveline files, elmer result files, or pvtu files
- `elmerrunlib.py` - run Elmer solver file 
- `floatlib.py` - calculate height above flotation
- `flowparameterlib.py` - calculate Glen's flow law parameter from a given ice temperature
- `fluxlib.py` - calculate ice flux through a fluxgate
- `geodatlib.py` - read velocities from [NASA MEaSUREs Greenland Ice Velocity from InSAR](https://nsidc.org/data/NSIDC-0481/versions/1)
- `geotifflib.py` - read and write geotiff files
- `glaclib.py` - load glacier flowlines and extents for different dates
- `inverselib.py` - guess a sliding law coefficient using the shallow ice approximation, which is used as the initial condition for basal inversions
- `masklib.py` - create masks based on the ice extent
- `meshlib.py` - create input files for developing finite element meshes in Gmsh and ElmerGrid 
- `picdatalib.py` - digitize and pick data from a figure from another paper
- `shapefactorlib.py` - calculate a shapefactor to account for lateral drag when using a flowline in Elmer


## big3/modeling
- `elmerlib`: functions and modules developed for [Elmer](https://github.com/ElmerCSC/elmerfem)
- `inversions`: scripts to run, process, and plot inversions
- `meshing`: scripts for writing meshes and outputting files necessary to run simulations in Elmer
- `pleiades_scripts`: scripts for running Elmer on [NASA Pleiades](https://www.nas.nasa.gov/hecc/#url)
- `simulations`: scripts for running diagnostic and prognostic simulations
- `solverfiles`: all Elmer solverfiles

## big3/observations
- scripts for processing and plotting observations

## big3/presentations
- scripts for figures in oral and poster presentations at conferences


## Resulting Publications
[Kehrl, L.M., I. Joughin. D.E. Shean, D. Floricioiu, and L. Krieger. 2017. Seasonal and interannual variabilities in terminus position, glacier velocity, and surface elevation at Helheim and Kangerlussuaq Glaciers from 2008 to 2016. Journal of Geophysical Research: Earth Surfaces, 122, doi:10.1002/2016JF004133](http://onlinelibrary.wiley.com/doi/10.1002/2016JF004133/full)
- Figures 1a, 2a, 3, 4, 6, 7: [big3/observations/velocity/velocity_vs_terminus.py](https://github.com/kehrl/big3/blob/master/observations/velocity/velocity_vs_terminus.py)
- Figures 1bc, 2bc: [big3/observations/elevation/elevation_along_flowline.py](https://github.com/kehrl/big3/blob/master/observations/elevation/elevation_along_flowline.py)
- Figures 5, 8, 9: [big3/observations/velocity/plot_spatial_variability.py](https://github.com/kehrl/big3/blob/master/observations/velocity/plot_spatial_variability.py)
- Figure 10: [big3/observations/elevation/thinning_fluxgate_dems.py](https://github.com/kehrl/big3/blob/master/observations/elevation/thinning_fluxgate_dems.py)
