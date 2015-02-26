# Defaults

The usual Go coding conventions apply unless overridden hereunder.
Refer to [Effective Go](https://golang.org/doc/effective_go.html) for
some advice.

## File Names

As far as possible, and is reasonable, use single word names for
files.  Where it is necessary to have more than one word, the
following apply.

- For source code files, separate the words using an **underscore**.
  E.g. `atom_builder.go`.
- For documentation files, separate the words using a **hyphen** (this
  file is an example).  E.g. `molecule-normal-form.md`.

In all cases, file names should be only in **lowercase English**.
They should be in plain ASCII, for ease of use within the terminal, as
well as for easy portability across operating systems and locales.

## Type Names

Types which are unexported should still use an initial uppercase
letter, except that they are prefixed with an underscore.
E.g. `_Atom`.

## Variables

For both top-level variables and `struct` members, in general, the
last part of a variable's name should indicate *what* it is.  This
applies to function and method parameter names as well.

```{go}
var numSingleBonds int  // Incorrect; this is NOT a collection of bonds.
var singleBondCount int // Correct; this IS a count.
```

In addition, variables should be named also according to the following
rules.

### Unexported Variables

Choose brief names for unexported variables.  Use common abbreviations
where applicable.  However, make sure that each name is unambiguous
and not cryptic.  Also, the abbreviations should not morph the names
into tongue twisters.

### Exported Variables

Names of exported variables should as full as necessary.  While you
can use common abbreviations, avoid less known ones.  However, avoid
needlessly long and overly descriptive names.

### Local Variables

Do *not* use the `var` form for local variables; use `:=`, and let Go
do the type inference.

## Functions and Methods

Use top-level functions where they are the best fit.  Do not create
`struct`s only to convert such functions into methods.

Top-level functions and methods follow the same set of conventions.
Both exported and unexported functions and methods follow the same set
of conventions.

### Names

Functions and methods representing actions should have names that are
verb forms, preferably in the imperative.

```{go}
func ringDetection() error { } // Incorrect; uses noun form.
func detectRings() error { }   // Correct; uses imperative verb form.
```

Functions and methods whose return value(s) (not side effects) are
primary, should have their names as noun words representing the nature
of the return values.

```{go}
func (m *Molecule) findBondBetween(atomId1, atomId2 uint16) *Bond { } // Incorrect.
func (m *Molecule) bondBetween(atomId1, atomId2 uint16) *Bond { }     // Correct.

```

### Return Values

Pay attention to what each function or method answers to its caller.
Choose carefully if they should be returned by value or by reference.

When multiple values are returned, if one of them represents a boolean
status or an error status/value, put that at the end of the return
value tuple.  This is in line with the 'comma ok' convention.

When a value and an error are returned, the following rules apply.

- **(good value, `nil`)** : best case; caller should continue.
- **(good value, non-`nil`)** : minor issue to note, but caller should
  continue.
- **(invalid value or `nil`, non-`nil`)** : major issue; caller should
  handle the error.
- **(invalid value or `nil`, `nil`)** : forbidden combination.

### Chaining

Where chaining function or method calls makes natural sense, the names
of such functions or methods should be as fluent as possible.
