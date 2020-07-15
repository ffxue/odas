# About ODAS

The Optimization-based Detection of Architectural Symmetries (ODAS) is a free, fast, accurate, and robust library for the problem of architectural symmetry detection from 3D point clouds.

# How to cite

  Xue, F., Lu, W., Webster, C., and Chen, K. (2019). An optimization-based approach for detecting architectural symmetries from 3D point clouds. Submitted to ISPRS Journal of Photogrammetry and Remote Sensing, 128, 32-40. doi: [10.1016/j.isprsjprs.2018.12.005](//doi.org/10.1016/j.isprsjprs.2018.12.005)

# How does it work

Symmetry in 3D space is fully determined by three real parameters, e.g., the plane (rho, heading, tilt) for a reflective symmetry, i.e., the symmetry axis; the 3D line (a, b, c) for a 3D rotation symmetry. For buildings and infrastructure, very often the symmetry axes are either perfectly horizontal or perfectly vertical, where only two parameters are involved. Thus, symmetry detection can be achieved by optimizing the parameters for the maximum point correspondence or the minimum RMSD (root-mean-square distance) of the symmetric transformation. As its name, ODAS incorporates modern optimization algorithms to solve the "parameter optimization" problem to detect the symmetries for buildings and infrastructures. In ODAS, over 20 algorithms are accessible. The DIRECT, CMA-ES, MLSL-LDS, and PSO are among the best options.

For more details, please refer to the papers.

# How to use ODAS

If you have a large-scale 3D point cloud, please convert the format to .pcd.

# Goal of the research project

The project is [A derivative-free optimization (DFO) approach to architectural symmetry detection from 3D point clouds](//www.researchgate.net/project/A-derivative-free-optimization-DFO-approach-to-architectural-symmetry-detection-from-3D-point-clouds).

This research project aims to recognize architectural symmetry from 3D point clouds based on known architectural styles and building technologies. Examples of architectural styles are neoclassical architecture and modernist architecture, which have been reflected in real-life physical architecture and influenced by architecture preferences and construction standards and other factors. Examples of building technologies include popular materials and building codes at the moment of construction.

The key scientific question is "How can the symmetry be efficiently recognized with constraints of architectural styles and knowledge by computer program and then be effectively integrated into as-built model generation methods?"

# Dependencies

[nlopt](//nlopt.readthedocs.io/)

[popot](//github.com/jeremyfix/popot)

[libcmaes](//github.com/beniz/libcmaes)

[nsga2-cpp](//github.com/dojeda/nsga2-cpp)

[Eigen 3.2+](//eigen.tuxfamily.org)

# Install

* Ubuntu 16.04, GCC 4.6+ 
* Install the dependencies listed above

    cmake .
    
    make -j2

# How to contribute

# License

LGPL-3.0

# Acknowledgements

This work was supported by 
* The Hong Kong Research Grant Council grant numbers [17200218](//cerg1.ugc.edu.hk/cergprod/scrrm00542.jsp?proj_id=17200218) and [17201717](//cerg1.ugc.edu.hk/cergprod/scrrm00542.jsp?proj_id=17201717) and
* The University of Hong Kong grant numbers 201702159013 and 201711159016
 
