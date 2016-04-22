#HTSOHM-dev

This is a package which can be used to map the structure-property space   
spanned by hypothetical, porous <i>pseudo materials</i>. This package contains   
modules that manage randomly-generating libraries of pseudo-materials,   
distributing many simulations across a computing cluster, and selectively   
mutating rare materials--in an iterative process.
 
## Getting started

The HSTOHM-dev libraries were written predominantly in Python3 and Bash. It   
also requires `RASPA-2.0`, available :   

https://github.com/numat/RASPA2   

To configure HTSOHM, you may run `configure.sh`, which writes the following    
to `~/.bashrc`:    
```
 export RASPA_DIR=${HOME}/RASPA        # for a default RASPA build
 export HTSOHM_DIR=${HOME}/HTSOHM-dev  # if cloned into ${HOME}
 export RASPA_DIR=${RASPA_DIR}
 export SRC_DIR=${HTSOHM_DR}/bin
 export FF_DIR=${RASPA_DIR}/share/raspa/forcefield
 export MAT_DIR=${RASPA_DIR}/share/raspa/structures/cif
```
To run the program, execute:    
  `python $HTSOHM_DIR/bin/HTSOHM.py <MpG> <NoA> <S0> <NoB>`    
`<MpG>` : number of materials per generation.   
`<NoA>` : number of atom-types.   
`<S0>`  : initial strength-parameter. 
`<NoB>` : number of bins.

By default the program will terminate after 20 generations.   
   
Please send questions/comments/concerns to `ark111@pitt.edu`.

## License

HTSOHM-dev and related modules released under the MIT License 2016. 

## Development environment

```
3.5.1 |Anaconda 2.4.1 (64-bit)| (default, Dec  7 2015, 11:16:01) 
[GCC 4.4.7 20120313 (Red Hat 4.4.7-1)]

Bash 4.1.2(1)-release
RHEL 6.6	2.6.32-358
```

