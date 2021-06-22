# AutoMapper - Beta

AutoMapper is a package of tools designed to automate the creation of files required for the [**LAMMPS**](https://lammps.sandia.gov) command _fix bond/react_. These tools allow the user to convert LAMMPS 'read_data' format files into partial structure molecule files which can be used to create a map file of a reaction. AutoMapper has been designed to work alongside [**Moltemplate**](https://github.com/jewettaij/moltemplate) output files, but will work in general cases so long as some assumptions are followed. See the _Assumptions_ section for details.

## Tool Description
All tools are accessible through `AutoMapper.py`. They are:
- `clean`: Change the types of a series of files to the lowest possible values whilst ensuring all atom types are included in each file. This is designed to work with Moltemplate which generates a list of all possible coefficients in a `systems.in.settings` file.
- `molecule`: Convert a LAMMPS `read_data` input file to the LAMMPS molecule file format.
- `molecule-partial`: Same as `molecule` except the output molecule file is reduced to atoms that are 3 or fewer bonds away from the bonding atoms. Edge atoms are found automatically and recorded in the file.
- `lammps-partial`: Same as `molecule-partial` except the output file is in the LAMMPS `read_data` format. Used to check partial structures are correct with visualiser tools like Ovito.
- `map`: Create a map from a pre-bond molecule file and a post-bond molecule file.

## General Usage

`AutoMapper.py` can be called from the command line with arguments for user control. It is recommended to put this downloaded repository on your `PATH` so that `AutoMapper.py` can be called from your modelling directories. Use `AutoMapper.py -h` for a full description of the arguments. Not obvious arguments are:
- `--ba: The atomIDs of the bonding atoms`
- `--da: The atomIDs of the atoms to be deleted, if any`
- `--ebt: A list of element symbols by atom type`

Below are command line examples for the different tools.

```
AutoMapper.py . clean pre-reaction.data post-reaction.data --coeff_file system.in.settings
AutoMapper.py . molecule cleanedpre-reaction.data --save_name pre- --ba 1 4
AutoMapper.py . molecule-partial cleanedpost-reaction.data --save_name pre- --ba 3 6 --da 7 9 --ebt H C C O
AutoMapper.py . lammps-partial cleanedpost-reaction.data --save_name post- --ba 3 6 --ebt H Si O O
AutoMapper.py . map pre-molecule.data post-molecule.data --ebt H H C C N O O 
```

## Typical Workflow
Assuming you have used Moltemplate, the inital files you will need are a LAMMPS file that represents the molecule structure before bonding (the 'pre-bond' file) and one representing the molecule structure after bonding (the 'post-bond' file). It is better that pre-bond and post-bond files only include the two molecules required for a reaction (or one molecule, in the case of an intramolecular reaction). You will also need the `system.in.settings` file from Moltemplate, or an equivalent file that includes all the required coefficients for LAMMPS. You can then use `clean` to remove unnecessary coefficients from the files, in a similar fashion to `cleanup-moletemplate.sh` from Moltemplate.

Then run `molecule` or `molecule-partial` on both new LAMMPS format files which will have the prefix 'cleaned'. For small molecules you may use `molecule` but it is recommended by fix bond/react to use the smallest partial molecule structure possible, which `molecule-partial` should find. If you wish to check the partial structure you may use the `lammps-partial` tool with the same arguments as `molecule-partial` to create a LAMMPS `read_data` format partial structure file. You may find this easier to load into visualiser tools like Ovito. Apart from save names, the only information the user must supply are the bonding atoms used in the reaction, and the element symbols by atom type in the case of the partial structure tools. The order of the bonding atoms is important: the first ID given must represent the same atom in the pre-bond and post-bond files, even if the ID number is different. The atoms to delete should be specified here if there are any atoms to be deleted in the reaction. 

Finally, the `map` tool can be used to create a map between the pre-bond and post-bond molecule files created previously. The only additional information needed is the element symbols by atom type series; all other information is taken from comments at the top of the molecule files. This will create a map file with all sections filled out. Your terminal may tell you to check some atom pairs if the pair has been found by inference. If you want to see how the path search has determined each mapped pair, use the `--debug` argument to print this information to your terminal.

## Assumptions
These tools have been built for the LAMMPS atom type 'full' and will **not** work with other atom types.
The `clean` tool requires all files have the same unit cell dimensions and that pair coefficients are defined with numbers and not wildcards (*)
The header of the LAMMPS data files must look as below for this code to read the required information. This is the standard output from Moltemplate. This assumption may be removed in the future.
```
# Comment - Either one line or multiple lines with a '#' preceding

a  atoms
b  bonds
c  angles
d  dihedrals
e  impropers

f  atom types
g  bond types
h  angle types
i  dihedral types
j  improper types

x1 x2 xlo xhi
y1 y2 ylo yhi
z1 z2 zlo zhi
```

## Limitations
Any byproducts created by the reaction (e.g. a water molecule) must be specified as delete atoms with the `--da` argument if using the `-partial` tools. Use the `molecule` tool if you don't wish to do this.
Further testing of more complicated molecules to come.