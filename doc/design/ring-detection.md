# Detection of Rings and Ring Systems

Accurate detection of a molecule's rings and ring systems is a
surprisingly complex problem.  Literature shows that it has been a
much debated and contentious subject as well!

In **RxnWeaver**, we attempt at striking a balance between theoretical
accuracy (which is debated) and practical accuracy (which determines
which reactions can take place and which cannot).

The following is a high-level description of the algorithm we employ.

## Ring Detection

1. Calculate the `frerejacque` number: `no. of bonds - no. of atoms + 1`.
1. If `frerejacque <= 0`, the molecule is acyclic.  Exit.
1. For now, we reject molecules such as fullerenes that have too many
   rings.  This may change in future.
1. Make a copy of the adjacency lists, _etc_.
1. Prune all terminal atom chains, so that we are left with only
   cyclic atoms.
1. At this point, if all the remaining atoms have exactly two
   neighbours each, there is only one ring.  Note that as the only
   ring, and as the only ring system.  Exit.

If we have come this far, there are at least two rings in the
molecule.  We now perform a breadth-first search, starting with a
non-junction atom.

1. Initialise a list of candidate paths.
1. Locate a random non-junction atom.  If no such exists, pick any
   atom.
1. Create a new path, and add this atom to it.
1. Add this path to the list of candidate paths.
1. While this list is non-empty, evaluate each path.

### Path Evaluation

1. For each neighbour of the most recent atom visited, the following
   is followed.
1. If the next atom is the same as the beginning of the path, then we
   have a candidate ring.  We validate it, and - if successful - add
   it as a ring to the molecule.  We proceed to the next neighbour.
1. If the next atom is one of the atoms visited already, then we have
   a candidate ring.  We validate it, and - if successful - add it as
   a ring to the molecule.  We proceed to the next neighbour.
1. Else, construct a new candidate path by appending the next atom to
   the current path.
1. Add the new candidate path to the list of candidate paths.

### Ring Validation

If the candidate ring has only three members, it is always a valid
ring (there is exactly one way of forming such a ring).

Otherwise, we employ the following procedure.

1. For each atom in the candidate ring, count the number of its
   neighbours that are present in this candidate ring.
1. If that number is exactly two, this ring is genuine.
1. If it is greater than two, _i.e.,_ at least one atom has more than
   two neighbours in this ring, then we have a spurious outer shell
   ring, not a genuine ring.

We evaluate all candidate rings, resulting in the set of genuine rings
in the molecule.

### Ring System Detection

We then detect ring systems formed by the rings.

- If a ring and another share a bond, they are fused.  Create a ring
  system of them.
- Each other ring sharing a bond or an atom with this ring system
  joins it.
- Disjoint rings start independent ring systems.
- Ring systems that share a ring are fused into a single ring system.
- Loop through the rings, until all of them fall into their respective
  ring systems.

Note that it is possible for a ring system to be monocyclic.
