# WiggleToSite
Adjusts glycosidic linkages to bring residues into a (second) binding site, to check if bidendate binding is possible.
Still very much scientist software at this point. May become a webtool at glycam.org.

# Dependencies:
GMML (https://github.com/GLYCAM-Web/gmml.git). Probably the dev version.

# Compile:
Run ./compile.sh

# Usage:
Run the executable: bin/WiggleToSite for instructions. See tests/ folder for example inputs.

# Input file:
I plan to have a nice user interface that writes this file. For now, it's difficult to understand. 
Here is an example:
https://github.com/gitoliver/WiggleToSite/blob/master/tests/input.txt

# Input file oddities:
The END lines are essential. 

The section naming must match pefectly (e.g. always "MovingResidues", never "MovingResidue")

Residues in your PDB files may have a chain ID (A,B,X etc). If they don't use a "?". 

Order the Moving and Target Residues in order that they should superimpose. Make sure the residue names and atom names match exactly (This may go away in furture versions).

List the linkages in the order you want them adjusted initially (it becomes random after the first iteration).

Specify each linkage from non-reducing to reducing (e.g. Gal1-4Glc, not Glc4-1Gal). The default is to move the atoms on the non-reducing (Gal) side of the example. To change that, add the "noReverse" keyword.

I know it's odd that it's not automagically figured out for you which residues can move, but this level of control allows you to move linkages in separate glycans on the same glycoprotein.
