# Threading

The construction, analysis and transformations of any given molecule
should be within a single thread whenever possible.

The following scenarios describe the notable exceptions.

# Bulk-processing

As the input representation of a molecule is read, a factory
constructs the in-memory molecule.  When bulk-processing molecules,
the input text is read by a dedicated reader thread, while batches of
input molecules are processed by separate worker threads.

# Retro-synthesis

In the case of a single goal molecule being tried for retro-synthesis,
each candidate reaction can be tried in a new thread.  Since the above
can apply to an arbitrary intermediate molecule too, for efficiency, a
pool of workers handles the generation of reactants for product
molecules.

In order to facilitate the above, a queue of product molecules is
maintained by a master thread.  Whenever a worker thread becomes
available, it picks up the next molecule to be processed off the said
queue.
