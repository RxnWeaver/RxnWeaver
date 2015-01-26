# General

**_N.B._ In what follows, we always assume a right-handed cartesian
  coordinate system.**

The following order of precedence is used in utilising the available
information in the input molecule, while determining stereo parities.

1. If 3-D coordinates are available for all molecules, or for all
   molecules except one, we use them for all calculations.
1. If 2-D coordinates are available for all molecules, and 2-D bond
   stereo information is available, we adjust Z-coordinates of the
   involved atoms appropriately.  Then, we employ the same
   calculations as above.
1. Should the stereo configuration involve only planar elements
   (stereogenic bond case), we utilise the 2-D coordinates of the
   atoms if available.

**_N.B._** Should the input molecule have *no* recorded stereo
information at all (_i.e._ all coordinates are zeroes or all relevant
Z-coordinates happen to be the same, up/down bond information is not
present, and identification of all chiral centres is not there), we
treat the molecule as having **no** stereo configuration at all!  The
reason is: in a general case, it is not possible to determine stereo
configuration unambiguously based only on topology and either 0-D or
2-D coordinates.

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

In this configuration, there is no specific atom acting as a stereo
centre.  Rather, the bond between **A** and **B** is assigned the
determined stereo parity.

### Parity Computation

Given the above configurations, and assuming **A** to have higher
priority than **B**, the following matrix determinant is calculated.

```
          | 1.0  X_x  X_y  X_z |
          | 1.0  A_x  A_y  A_z |
          | 1.0  B_x  B_y  B_z |
          | 1.0  Y_x  Y_y  Y_z |
```

Should the determinant be positive, the parity is `EVEN`; should it be
negative, it is `ODD`.

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

Again, in this configuration, there is no specific atom acting as a
stereo centre.  Rather, the bond between **A** and **B** is assigned
the determined stereo parity.

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

**_N.B._** *Should a chiral centre have two neighbours that are both
tetrahedral chiral centres themselves, the one with `EVEN` parity is
given priority over that with `ODD`.*

### Adding Z-coordinates

Next, we look at each tetrahedral chiral centre.  For illustration,
let us consider one such, which we refer to as **X**.  We designate
its neighbours **A**, **B**, **C** and **D**, in **_descending_**
order of priority.

Where 3-D coordinates are not available, we assume the input to
specify 2-D coordinates of each atom, with the Z-coordinate being
`0.0`.  For each of its bonds, the input is expected to specify
whether it is *on the plane*, *out of the plane* or *into the plane*.
The pointy end of the bond is always directed towards **X**.

Now, we walk through the bonds that are either out of the plane or
into the plane.  The end-point neighbour of such a bond is adjusted to
have its Z-coordinate either incremented by or decremented by a
pre-computed amount, respectively.

### Parity Computation

Given the above orders of priority of the neighbours **A**, **B**,
**C** and **D**, we calculate the following matrix determinant.

```
          | 1.0  A_x  A_y  A_z |
          | 1.0  B_x  B_y  B_z |
          | 1.0  C_x  C_y  C_z |
          | 1.0  D_x  D_y  D_z |
```

Should the determinant be positive, the parity is `EVEN`; should it be
negative, it is `ODD`.

### Tetrahedral Atom With 3 Neighbours

In the case of a central tetrahedral atom with either one double bond
or an implicit hydrogen, we construct the matrix with the central atom
as that with the highest priority, followed by the three neighbours in
decreasing order of priority.

Given **X** as the central atom, and **A**, **B** and **C** as its
three neighbours listed in descending order of priority, we compute
the following matrix determinant.

```
          | 1.0  X_x  X_y  X_z |
          | 1.0  A_x  A_y  A_z |
          | 1.0  B_x  B_y  B_z |
          | 1.0  C_x  C_y  C_z |
```

Should the determinant be positive, the parity is `EVEN`; should it be
negative, it is `ODD`.

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
priority than **B**, the following matrix determinant is calculated.

```
          | 1.0  X_x  X_y  X_z |
          | 1.0  A_x  A_y  A_z |
          | 1.0  B_x  B_y  B_z |
          | 1.0  Y_x  Y_y  Y_z |
```

Should the determinant be positive, the parity is `EVEN`; should it be
negative, it is `ODD`.
