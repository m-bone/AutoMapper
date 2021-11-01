# AutoMapper

AutoMapper is a package of tools designed to automate the creation of files required for the [**LAMMPS**](https://lammps.sandia.gov) command _fix bond/react_. The tools allow the user to convert LAMMPS input files into molecule format files and an automatically generated a map file of a reaction. AutoMapper will reduce these files to the smallest possible partial structure without any additional user input. AutoMapper should work with any valid LAMMPS input file, though was designed to work with [**Moltemplate**](https://github.com/jewettaij/moltemplate) output files so users may find it easier to use these packages together. AutoMapper is completely forcefield independant so users are free to use any forcefield they wish. Any problems encountered while using AutoMapper should be raised in the [**Issues**](https://github.com/m-bone/AutoMapper/issues) section where they will be addressed as quickly as possible. For a full tutorial and annotated walkthrough, see the [**manual**](https://github.com/m-bone/AutoMapper/blob/Dev-Branch/AutoMapper_Manual.pdf).

Note: AutoMapper was built for LAMMPS version 29Oct20. The REACT package has been changed in newer versions, but it is believed AutoMapper still works. This will be fully tested in the future.

## Tool Description
All tools are accessible through `AutoMapper.py`. They are:
- `clean`: Unify the types (e.g. Atom, Bond, Angle, etc.) between two or more files and remove unused coefficients.
- `molecule`: Convert a LAMMPS input file to the LAMMPS molecule file format.
- `map`: Create pre- and post-bond molecule files and a map file, the full requirements for using `bond/react`.

## General Usage

`AutoMapper.py` can be called from the command line with arguments for user control. It is recommended to put this downloaded repository on your `PATH` so that `AutoMapper.py` can be called from your modelling directories. Use `AutoMapper.py -h` for a full description of the arguments or see the [**manual**](https://github.com/m-bone/AutoMapper/blob/Dev-Branch/AutoMapper_Manual.pdf) for more details. A full [**Example Workflow**](https://github.com/m-bone/AutoMapper/tree/main/Example_Workflow) is available within the repository which walks through the use of AutoMapper within a Linux environment.

Below are command line examples for the different tools.

```
AutoMapper.py . clean pre-reaction.data post-reaction.data --coeff_file system.in.settings
AutoMapper.py . molecule cleanedpre-reaction.data --save_name pre- --ba 1 4
AutoMapper.py . map cleanedpre-reaction.data cleanedpost-reaction.data --save_name pre-molecule.data post-molecule.data --ba 2 5 3 7 --da 1 8 1 8 --ebt H H C C N O O 
```

## Assumptions
AutoMapper requires Python 3.6+ and the third party Python module [**natsort**](https://pypi.org/project/natsort/) to run.
These tools have been built for the LAMMPS atom style 'full'; results with other atom styles may vary. 
For the `clean` tool pair coefficients must be defined with numbers and not use wildcards (`*`).
