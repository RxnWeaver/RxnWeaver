# Aromaticity Determination

Determining whether a ring or a ring system is aromatic follows the
following rules.

- Aromaticity of a ring system is given priority over that of its
  constituent rings.
- Accordingly, when a ring system, as a whole, is aromatic, the
  individual rings in it are *not* specifically marked as aromatic.
- On the other hand, when a ring system is *not* aromatic, we
  investigate each constituent ring individually.  Each such ring
  found to be aromatic is marked as such.
- In all cases, the atoms and bonds of an aromatic system are all
  marked as being part of one such.

When the entirety of a molecule is covered by exactly one ring system,
and that ring system is aromatic, we say that the molecule itself is
aromatic.

## Hückel's Rule

Application of Hückel's rule requires the computation of the number of
total number of delocalised π-electrons in the relevant ring system or
ring.  In order to facilitate that, each atom keeps track of its
current π-electron count.

The bonds an atom participates in, and the residual charge on the
atom - if any - has to be duly taken into account for a correct
calculation of delocalised π-electrons.

Hückel's rule is a convenient way of determining if a planar ring
system or ring is aromatic.  The rule says that the number of
delocalised π-electrons in the ring system or ring should be even, but
not a multiple of `4`.  It can, therefore, be represented by the
following formula.

```
        no. of π-electrons = 4 * n + 2, n = 0, 1, 2, 3, 4, 5, 6
```

The correctness of Hückel's rule does not appear to hold beyond that
limit.  Even so, several exceptions to the rule are known.  They need
special handling.
