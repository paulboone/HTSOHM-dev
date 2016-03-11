#HTSOHM-dev

This is a package which can be used to map the structure-property space   
spanned by hypothetical, porous <i>pseudo materials</i>. This package contains   
modules that manage randomly-generating libraries of pseudo-materials,   
distributing many simulations across a computing cluster, and selectively   
mutating rare materials--in an iterative process.
 
## Getting started

The HSTOHM-dev libraries were written predominantly in Python3 and Bash. It   
also requires `RASPA-2.0`, available : <find Dubbledam link>.   
In order to run `workflow-*.ipynb` one must also have the Jupyter Notebook   
installed (`$ pip install notebook`).   
   
In order to run this package as-is, the users `$HOME` directory must be   
contain the following:   
```
 ${HOME}/RASPA/bin/simulate
 ${HOME}/RASPA/share/raspa/structures/cif/
 ${HOME}/RASPA/share/raspa/forcefield/
 #created by default RASPA installation

 ${HOME}/HTSOHM-dev/data/
 ${HOME}/HTSOHM-dev/bin/*
 ${HOME}/HTSOHM-dev/workflow-*.ipynb
 #after git cloning repo into $HOME 
```
   
A Jupyter notebook `workflow-*.ipynb` walks the user through the various   
steps involved in this computational method. To get started:   
   
`$ jupyter notebook workflow-*.ipynb`
   
Please send questions/comments/concerns to `ark111@pitt.edu`.

## License

HTSOHM-dev and related modules released under the GNU General Public License v3.0. 

## Development environment

```
3.5.1 |Anaconda 2.4.1 (64-bit)| (default, Dec  7 2015, 11:16:01) 
[GCC 4.4.7 20120313 (Red Hat 4.4.7-1)]

Bash 4.1.2(1)-release
RHEL 6.6	2.6.32-358
```

