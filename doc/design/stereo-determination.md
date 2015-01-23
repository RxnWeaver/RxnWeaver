# Stereo Configuration Determination

The major types of stereo configurations have to be determined when
analysing the structure of a molecule.

- Stereogenic bond parity
- Cumulene stereo parity
- Tetrahedral stereo parity
- Allene stereo parity
- Alternating stereo bonds in rings

Our treatment of the same shall follow the approach described in the
InChI API documentation v1.04.  Specific algorithms may coincide or
differ, based on computational considerations.

Notably, we differ in parity naming convention.  InChI maps all
determinable configurations into two parities: **even** and **odd**.
This is not very intuitive, and often requires looking up the meaning
of a parity in the context of a particular stereo configuration.  We
have observed that most scientists re-map the InChI parity value to
_old school_ `R`, `S`, `E`, `Z`, _etc_.  Accordingly, we directly use
the old school symbols themselves.

**_N.B._** Should the input molecule have *no* recorded stereo
information at all (_i.e._ up/down bond information and identification
of all chiral centres), we treat the molecule as non-chiral!  The
reason is: in a general case, it is not possible to determine stereo
configuration unambiguously based only on topology and 2-D
coordinates.

## Bond Vector Construction

Let **O** be the origin of the coordinate system.  The 3-D coordinates
`(x, y, z)` of an atom **P** represent a 3-D vector **_OP_**.  For
illustration, suppose that a bond exists between atoms **A** and
**B**.  We can construct the bond vector between them **_AB_** as

- **_AB_** = **_OB_** **-** **_OA_**.

## Stereogenic Bond Parity

Let atoms **A** and **B** be singly or doubly bonded.  Suppose that
atom **X** is singly bonded to **A**, and that atom **Y** is singly
bonded to **B**.  Further suppose that all the bonds lie in the same
plane.  There are two possible configurations, as shown here.

```
    X                      X      Y
     \                      \    /
      A==B                   A==B    (or)    A==B
          \                                 /    \
           Y                               X      Y
```

We compute the vector dot product of **_AX_** and **_BY_**.

- **dprod** = **_AX_** **.** **_BY_**

Should **dprod** be negative, we have the former case with those two
vectors pointing in opposite directions, falling on opposite sides of
the stereogenic bond between **A** and **B**.  Such a configuration is
named `E` parity (for German *entgegen*, meaning opposite).

Should **dprod** be positive, we have the one of the latter cases with
those two vectors falling on the same side of the stereogenic bond
between **A** and **B**.  Such a configuration is named `Z` parity
(for German *zusammen*, meaning same).

## Cumulene Stereo Parity

The treatment of cumulenes follows that of stereogenic bonds described
above.  In the case of cumulenes, atoms **A** and **B** are connected
by a chain of double bonds.  The two possible stereo configurations
are as follows:

```
    X                           X           Y
     \                           \         /
      A==...==B                   A==...==B    (or)    A==...==B
               \                                      /         \
                Y                                    X           Y
```

where, `...` represents one or more atoms all of which are doubly
bonded to form a chain.

## Tetrahedral Stereo Parity

Here is an outline of how we determine tetrahedral stereo parity.

### Priority Determination

First, we determine the base priorities of all atoms.  The priority of
an atom is the sum of its atomic number and those of all of its
neighbours.  For a double bond, the neighbour is counted twice, _etc_.

An indicative function to perform the same could look as follows, in
Rust (here we assume that this method is in `Molecule`).

```{rust}
fn assign_priorities(&mut self: Self) {
    for &a in self.atoms() {
        a.set_priority(
            a.bonds()
            .fold(a.atomic_number(),
                  |acc, &bnd| { acc +
                                bnd.order() *
                                bnd.other_atom(a).atomic_number() })
        );
    }
}
```

*Should a chiral centre have two neighbours that are both tetrahedral
chiral centres themselves, the one with `R` parity is given priority
over that with `S`.*

### Adding Z-coordinates

Next, we look at each tetrahedral chiral centre.  For illustration,
let us consider one such, which we refer to as **X**.  We designate
its neighbours **A**, **B**, **C** and **D**, in **_descending_**
order of priority.  For each of its bonds, the input is expected to
specify whether it is *on the plane*, *out of the plane* or *into the
plane*.  The pointy end of the bond is always directed towards **X**.

We assume the input to specify the 3-D coordinates of each atom, where
the Z-coordinate is `0.0`.

Now, we walk through the bonds that are either out of the plane or
into the plane.  The end-point neighbour of such a bond is adjusted to
have its Z-coordinate either incremented by or decremented by `0.5`,
respectively.

### Stereo Parity Determination

An important matter to realise (and visualise) is that since **A**,
**B** and **C** lie in a plane, **D** falls to one side of that plane.

We compute the vector cross-product of **_XA_** and **_XB_**.

- **_XA_x_XB_** = **_XA_** **x** **_XB_**

Depending on the angle `theta` between them (as measured from
**_XA_**), the resulting vector **_XA_x_XB_** points towards either
that side of the plane with **D** in it, or the opposite side.

The case is decided by doing a dot product of this result vector with
**_XD_**.

- **dprod** = **_XA_x_XB_** **.** **_XD_**

Should **dprod** be positive, then **_XA_x_XB_** was a result of a
clockwise rotation of priority from **A** to **B**.  Else, it was
anti-clockwise.

As per convention, clockwise rotation parity is named `R` (for Latin
*rectus*, meaning right) and anti-clockwise rotation parity is named
`S` (for Latin *sinister*, meaning left).

### Simplified Computation

Computationally, the above can be reduced to the determinant of the
following 3x3 matrix (or its row-major equivalent).

```
| XA_x  XB_x  XD_x |
|                  |
| XA_y  XB_y  XD_y |
|                  |
| XA_z  XB_z  XD_z |
```

The equivalence of the two methods of calculations follows directly
from the expansion of the vector computations, and mapping of the
terms to those in the calculation of determinants using Sarrus' rule.

Should the value of this determinant be positive, the parity is `R`;
else it is `S`.
