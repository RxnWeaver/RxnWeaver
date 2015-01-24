# Stereo Configuration Determination

The major types of stereo configurations have to be determined when
analysing the structure of a molecule.

- Stereogenic bond parity
- Cumulene stereo parity
- Stereo bonds in rings
- Tetrahedral stereo parity
- Allene stereo parity

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

where, `...` represents two or more atoms all of which are doubly
bonded to form a chain.

## Stereo Bonds in Rings

As per InChI rules, in all-carbon rings whose size >= 8, all double
bonds are treated as stereogenic.

Similarly, in all-carbon rings whose size >= 8, all alternating
single/double bonds are treated as stereogenic.  That helps in
unambiguously distinguishing between the following two ring
configurations.

```
                 D--E
                //   \\
             B--C      F--G                 B--C      F--G
           //              \\             //   \\   //    \\
          A                  H           A      D--E        H
           \               /              \               /
            N==M       J==I                N==M       J==I
                \     /                        \     /
                  L==K                           L==K

```

For each stereogenic bond, the determination of stereo parity shall be
according to the same rule given above for simple stereogenic bonds.

## Tetrahedral Stereo Parity

Here is an outline of how we determine tetrahedral stereo parity.

### Priority Determination

First, we determine the base priorities of all atoms.  The priority of
an atom is the sum of its atomic number and those of all of its
neighbours.  For a double bond, the neighbour is counted twice, _etc_.

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

### Tetrahedral Atom With 3 Neighbours

In the case of a central tetrahedral atom with either one double bond
or an implicit hydrogen, the following rule shall apply.

Suppose that atoms **A**, **B** and **C** are the neighbours of the
central atom **X**, in **_descending_** order of priority.

Now, **A**, **B** and **C** lie in a plane, with **X** falling to one
side of that plane.  The following calculations assume that **X** is
*behind* the said plane.

We compute the vector cross-product of **_XA_** and **_XB_**.

- **_XA_x_XB_** = **_XA_** **x** **_XB_**

Depending on the angle `theta` between them (as measured from
**_XA_**), the resulting vector **_XA_x_XB_** points towards either
that side of the plane with **X** in it, or the opposite side.

The case is decided by doing a dot product of this result vector with
**_XC_**.

- **dprod** = **_XA_x_XB_** **.** **_XC_**

Should **dprod** be negative, then **_XA_x_XB_** was a result of a
clockwise rotation of priority from **A** to **B**.  Else, it was
anti-clockwise.

As per convention, clockwise rotation parity is named `R` and
anti-clockwise rotation parity is named `S`.

## Allene Stereo Parity

The case of allene stereo configurations is similar to that of
cumulene ones.  The key difference is that in the case of cumulene,
there is no atom identified as a chiral centre.  The cumulene chain
represents the stereogenic bond.  On the other hand, in the case of
allene, the central `C` is marked as a chiral centre.

```
            X         Y
             \       /
              A==C==B
```

Here, **A** should have a higher priority than **B**.

Given the above arrangement, there are two possible stereo
configurations.  They are both determined based on how **X** and **Y**
are oriented relative to each the other, when observed from **A**,
with the line containing **A**, **C** and **B** collapsed.

```
                     Y                      Y
                     |                      |
                  X--A                      A--X
```

Looking at this in another way, given that the atom `C` is a chiral
centre, we can test if (**A**, **B**, **Y**) represents a clockwise
rotation or an anti-clockwise, when seen from **X** along **A**.  In
that case, the above two configurations translate to the following
two, respectively.

```
               Y                            Y
               |                            |
               B==C==A                A==C==B
```

In the first of the above cases, the parity is `EVEN`; it is `ODD` in
the latter case.

### Parity Computation

Given the above configurations, and assuming **A** to have higher
priority than **B**, the following matrix determinant should be
calculated.

```
          | 1.0  X_x  X_y  X_z |
          | 1.0  A_x  A_y  A_z |
          | 1.0  B_x  B_y  B_z |
          | 1.0  Y_x  Y_y  Y_z |
```

Should the determinant be positive, the parity is `EVEN`; should it be
negative, it is `ODD`.
