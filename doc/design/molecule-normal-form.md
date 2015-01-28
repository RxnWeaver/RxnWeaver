# Representation of Molecules

The internal representation format of molecules incorporates
information from both the input MDL `.MOL` V2000 format data, and the
additional information that **RxnWeaver** needs and updates throughout
the retrosynthesis process.

**_N.B._** _We do **not** provide compatibility with MDL's V3000
  format (`.MOL` or any other) files._

## Atoms

Each atom knows its element, current valence configuration, if it is
an isotope of its element, any residual charge on it, the number of
(explicit and implicit) hydrogen atoms attached to it, its spatial
coordinates, its neighbours, the specific bonds with them, _etc_.

Atoms have a normalised form to which they are converted after the
entire molecule is read and validated.  The process of normalisation
involves multiple steps, the most important of which are described
hereunder.

Each atom is assigned a priority as per the standard
**Cahn-Ingold-Prelog (CIP)** rules, authoritatively described in
*R.S. Cahn, C.K. Ingold and V. Prelog, Angew. Chem. 78, 413-447
(1966), Angew. Chem. Internat. Ed. Eng. 5, 385-415, 511 (1966); and
V. Prelog and G. Helmchen, Angew. Chem. 94, 614-631 (1982),
Angew. Chem. Internat. Ed. Eng. 21, 567-583 (1982)*.

Information about neighbours is stored as an adjacency list.  The
neighbours in the list are sorted in **descending** order of their
respective CIP priorities.

For ease of look up, adjacency lists are full: if atom **A** is a
neighbour of atom **B**, the adjacency list of **B** lists **A**, and
_vice versa_.  However, due to the CIP priorities mentioned, their
respective indices may be different.  This repetition does waste some
space, but is very convenient when accessing and traversing
neighbours.

## Bonds

Each bond represents the chemical bond between exactly one pair of
atoms.  It records the bond order, as well as stereo information
should it be applicable.

In a given molecule, we should never have more than one bond
connecting the same pair of atoms.

*N.B. We do not use InChI's provision to store stereo information
 relative to the second atom.  We always store stereo information of a
 bond relative to the first atom.*

In addition, each bond tracks if it is aromatic, if it participates in
any rings, which rings it participates in if so, if it is in a chain
linking two separate ring systems (that belong to the same connected
component), _etc_.
