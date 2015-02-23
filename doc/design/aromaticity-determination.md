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

## Some Typical and Atypical Examples
 References:
 Advanced Organic Chemistry, Reactions, Mechanisms and Structures, 4th Ed., Jerry March (AOC).
 
| Ring Name(Charge)          | Reference  |Pi Electrons|is Ring Aromatic?|is Ring System Aromatic? |
| :--------------------------|:-----------|:----------:|:---------------:|:-----------------------:|
| Triptycene                 |            |18          |PARTLY Y         |N                        |
| 18-annulene                |            |18          |Y                |Y                        |
| 14-annulene                |            |14          |Y                |Y                        |
| Anthracene                 | AOC43      |14          |N                |Y                        |
| Phenathrene                | AOC43      |14          |PARTLY Y         |Y                        |
| Phenalene                  | AOC44      |12          |PARTLY Y         |N                        |
| Phenalenide(-1)            | AOC44      |14          |N                |Y                        |
| Azulene                    | AOC49      |10          |N                |Y                        |
| Napthalene                 | AOC43      |10          |N                |Y                        |
| Benzene                    | AOC43      |6           |Y                |Y                        |
| 4H-pyran                   | AOC42      |6           |N                |N                        |
| Pyridine                   | AOC42      |6           |Y                |Y                        |
| Pyrylium(+1)               | AOC42      |6           |Y                |Y                        |
| 3‐oxathiazine‐2,2,4‐trione |            |6           |N                |N                        |
| 1,3-Cyclopentadienide(-1)  | AOC46      |6           |Y                |Y                        |
| 1,3-Cylopentadiene         | AOC46      |4           |N                |N                        |
| Furan                      | AOC45      |6           |Y                |Y                        |
| Pyrazole                   |            |6           |Y                |Y                        |
| Imidazole                  |            |6           |Y                |Y                        |
| Isothiazole                |            |6           |Y                |Y                        |
| Isoxazole                  |            |6           |Y                |Y                        |
| Oxazole                    |            |6           |Y                |Y                        |
| Thiazole                   |            |6           |Y                |Y                        |
| Thiophene                  | AOC45      |6           |Y                |Y                        |
| Pyrrole                    | AOC45      |6           |Y                |Y                        |
| Pyridine-N-oxide           | AOC42      |6           |Y                |Y                        |
| Pyridinium(+1)             | AOC42      |6           |Y                |Y                        |
    
