Stereo Configuration Determination
==================================

The major types of stereo configurations have to be determined when analysing the structure of a molecule.

- Tetrahedral chirality
- Double bond stereo parity
- Cumulene/allene stereo parity

Our treatment of the same shall follow the approach described in the InChI API documentation v1.04.  Specific algorithms may coincide or differ, based on computational considerations.

**N.B.** Should the input molecule have *no* recorded stereo information at all (*i.e.* up/down bond information and identification of all chiral centres), we treat the molecule as non-chiral!  The reason is: in a general case, it is not possible to determine stereo configuration unambiguously purely based on topology and 2-D coordinates.

Tetrahedral Chirality
---------------------

Here is an outline of how we determine tetrahedral chirality.

First, we determine the base priorities of all atoms.

```
assign_priorities_2(molecule) {
    molecule
    .atoms()
    .foreach((a) => {
        a.set_priority(
            a.neighbours()
            .map((nbr) => nbr.atomic_number())
            .sum()
            + a.atomic_number()
        );
    });
}
```

Next, we look at each tetrahedral chiral centre.  For illustration, let us consider one such, which we refer to as **X**.  We designate its neighbours **A**, **B**, **C** and **D**.  For each of its bonds, the input is expected to specify whether it is *on the plane*, *out of the plane* or *into the plane*.  The pointy end of the bond is always directed towards **X**.

Let **O** be the origin of the coordinate system.  We assume the input to specify the 3-D coordinates of each atom, where the Z-coordinate is `0.0`.  Thus, the 3-D coordinates `(x, y, z)` of an atom **P** represent a 3-D vector **_OP_**.