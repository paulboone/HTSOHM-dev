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

To configure HTSOHM, you may run `configure.sh`, which exports    
`$HTSOHM_DIR` to the directory in which the package was cloned.    
The RASPA directory (`$RASPA_DIR`) assumes a default installation   
in `$HOME`. The shell script also writes the following to `~/.bashrc`:    
```
 export HTSOHM_DIR=${HTSOHM_DIR}
 export RASPA_DIR=${RASPA_DIR}
 export SRC_DIR=${HTSOHM_DR}/bin
 export FF_DIR=${RASPA_DIR}/share/raspa/forcefield
 export MAT_DIR=${RASPA_DIR}/share/raspa/structures/cif
```
To run the program, execute:    
  `python $HTSOHM_DIR/bin/HTSOHM.py <MpG> <NoA> <S0> <NoB>`
Where `<MpG>`is the number of materials per generation, and    
`<NoA>` the number of atom-types, `<S0>` the initial strength-    
parameter, and `<NoB>` the number of bins.

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

