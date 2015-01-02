Synthesis Tree
==============

The retro-synthesis of each user-specified goal molecule constructs a synthesis tree.  The tree has the following properties.

- There exists exactly one node with no parent(s): *root* node, which represents the user-specified goal molecule.
- Every other node has at least one, possibly more, parent(s).
- Each node can have multiple children nodes.
- Leaf nodes represent those molecules that are commercially available, and hence do not need to be synthesised.
- Every node, other than leaf nodes, has one or more incoming product end-points.

End-points can be of two kinds:

- **linear synthesis end-point** : has one incoming path emanating from an upstream molecule, and
- **convergent synthesis end-point** : has two (or, rarely, more) incoming paths from as many distinct upstream molecules.

Consequently, it is possible for a given molecule to be the product of several linear and convergent syntheses.