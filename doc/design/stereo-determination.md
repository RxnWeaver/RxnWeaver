# Stereo Configuration Determination

The major types of stereo configurations have to be determined when analysing the structure of a molecule.

- Tetrahedral chirality
- Double bond stereo parity
- Cumulene/allene stereo parity

Our treatment of the same shall follow the approach described in the InChI API documentation v1.04.  Specific algorithms may coincide or differ, based on computational considerations.

**N.B.** Should the input molecule have *no* recorded stereo information at all (*i.e.* up/down bond information and identification of all chiral centres), we treat the molecule as non-chiral!  The reason is: in a general case, it is not possible to determine stereo configuration unambiguously purely based on topology and 2-D coordinates.

## Tetrahedral Chirality

Here is an outline of how we determine tetrahedral chirality.

### Priority Determination

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

### Adding Z-coordinates

Next, we look at each tetrahedral chiral centre.  For illustration, let us consider one such, which we refer to as **X**.  We designate its neighbours **A**, **B**, **C** and **D**, *in **descending** order of priority*.  For each of its bonds, the input is expected to specify whether it is *on the plane*, *out of the plane* or *into the plane*.  The pointy end of the bond is always directed towards **X**.

Let **O** be the origin of the coordinate system.  We assume the input to specify the 3-D coordinates of each atom, where the Z-coordinate is `0.0`.

Now, we walk through the bonds that are either out of the plane or into the plane.  The end-point neighbour of such a bond is adjusted to have its Z-coordinate either incremented by or decremented by 0.5, respectively.

### Bond Vector Construction

The 3-D coordinates `(x, y, z)` of an atom **P** represent a 3-D vector **_OP_**.  Thus, the current atoms of interest are represented by the vectors **_OX_**, **_OA_**, **_OB_**, **_OC_** and **_OD_**.  The vectors representing the bonds can then be constructed as follows.

- **_XA_** = **_OA_** **-** **_OX_**
- **_XB_** = **_OB_** **-** **_OX_**
- **_XC_** = **_OC_** **-** **_OX_**
- **_XD_** = **_OD_** **-** **_OX_**

### Stereo Parity Determination

An important matter to realise (and visualise) is that since **A**, **B** and **C** lie in a plane, **D** falls to one side of that plane.

We compute the vector cross-product of **_XA_** and **_XB_**.

- **_XA_x_XB_** = **_XA_** **x** **_XB_**

Depending on the angle `theta` between them (as measured from **_XA_**), the resulting vector **_XA_x_XB_** points towards either that side of the plane with **D** in it, or the opposite side.

The case is decided by doing a dot product of this result vector with **_XD_**.

- **dprod** = **_XA_x_XB_** **.** **_XD_**

Should **dprod** be positive, then **_XA_x_XB_** was a result of a clockwise rotation of priority from **A** to **B**.  Else, it was anti-clockwise.

As per convention, clockwise rotation parity is named `R` and anti-clockwise rotation parity is named `S`.
