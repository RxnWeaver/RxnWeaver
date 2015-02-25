package molecule

import (
	"fmt"
	"math"

	bits "github.com/willf/bitset"

	cmn "github.com/RxnWeaver/rxnweaver/common"
)

// Atom represents a chemical atom.
//
// It has various physical and chemical properties, apart from the
// transient information attached to it during its participation in
// reactions.
type _Atom struct {
	mol    *Molecule // Containing molecule of this atom.
	atNum  uint8     // Atomic number of this atom's element.
	symbol string    // Symbol, in case of a different isotope.
	iId    uint16    // Serial input ID of this atom.
	nId    uint16    // Normalised ID of this atom.

	X float32 // X-coordinate of this atom.
	Y float32 // Y-coordinate of this atom.
	Z float32 // Z-coordinate of this atom.

	hCount  uint8       // Number of implicit + explicit H atoms attached to this atom.
	charge  int8        // Residual net charge of this atom.
	valence int8        // Current valence configuration of this atom.
	radical cmn.Radical // Current radical configuration.

	unsaturation cmn.Unsaturation // Current composite state of this atom.

	pHash uint64 // A pseudo-hash of this atom, using some attributes.
	sHash uint64 // A pseudo-hash of this atom, using some attributes.

	bonds           *bits.BitSet // Bitmap of bonds of this atom.
	nbrs            []uint16     // Expanded list of neighbours of this atom.
	singleBondCount uint8        // Number of single bonds this atom has.
	doubleBondCount uint8        // Number of double bonds this atom has.
	tripleBondCount uint8        // Number of triple bonds this atom has.

	rings *bits.BitSet // Bitmap of IDs of rings this atom participates in.
	// Does this atom participate in at least one aromatic ring?
	isInAroRing bool
	// Is this atom a bridgehead of a bicyclic system of rings?
	isBridgeHead bool
	// Is this atom the sole common atom of all of its rings?
	isSpiro bool

	// The functional groups substituted on this atom.  They are listed in
	// descending order of importance.  The first is the primary feature.
	features []uint16

	// Number of electron-donating neighbours.
	edNbrCount int
	// Number of unsaturated electron-withdrawing neighbours.
	unsatEwNbrCount int
	// Number of saturated electron-withdrawing neighbours.
	satEwNbrCount int
}

// newAtom constructs and initialises a new atom of the given element
// type, and belonging to the given molecule.
func newAtom(mol *Molecule, atNum uint8, iId int) *_Atom {
	atom := new(_Atom)
	atom.mol = mol
	atom.atNum = atNum
	atom.iId = uint16(iId)

	el := cmn.PeriodicTable[cmn.ElementSymbols[atNum]]
	atom.symbol = el.Symbol
	atom.valence = el.Valence

	atom.bonds = bits.New(cmn.MaxBonds)
	atom.nbrs = make([]uint16, 0, cmn.MaxBonds)
	atom.rings = bits.New(cmn.MaxRings)

	atom.features = make([]uint16, 0, cmn.MaxFeatures)

	return atom
}

// AtomicNumber answers the atomic number of this atom.
func (a *_Atom) AtomicNumber() uint8 {
	return a.atNum
}

// Parent answers the parent molecule of this atom.
func (a *_Atom) Parent() *Molecule {
	return a.mol
}

// determineUnsaturation computes a composite metric that reflects the
// current state of the atom.
//
// This method is expected to be invoked during a molecule's
// normalisation only.
//
// Note that this method underlies - directly or indirectly - major
// parts of RxnWeaver's decision rules.  Exercise great caution should
// you need to modify this in any manner!
func (a *_Atom) determineUnsaturation() error {
	nb := int(a.bonds.Count())
	nn := len(a.nbrs)

	// Atom has a residual charge.
	if a.charge != 0 {
		a.unsaturation = cmn.UnsaturationCharged
		return nil
	}

	// For an uncharged atom, valence should be sane.
	if a.hCount > 0 {
		os := int8(nn) + int8(a.hCount)
		if ok, err := cmn.IsValidOxidationState(a.atNum, os); !ok {
			return err
		}
	}

	// Case of all single bonds.
	if nb == nn {
		a.unsaturation = cmn.UnsaturationNone
		return nil
	}

	// Double or triple bonds exist.
	ndb := 0
	nhdb := 0
	ntb := 0
	nhtb := 0
	mol := a.mol
	for bid, ok := a.bonds.NextSet(0); ok; bid, ok = a.bonds.NextSet(bid + 1) {
		b := mol.bondWithId(uint16(bid))
		oaid := b.otherAtomIid(a.iId)
		oa := mol.atomWithIid(oaid)
		switch b.bType {
		case cmn.BondTypeDouble:
			ndb++
			if oa.atNum != 6 {
				nhdb++
			}
		case cmn.BondTypeTriple:
			ntb++
			if oa.atNum != 6 {
				nhtb++
			}
		}
	}

	if ntb > 0 {
		if nhtb > 0 {
			a.unsaturation = cmn.UnsaturationTripleBondW
		} else {
			a.unsaturation = cmn.UnsaturationTripleBondC
		}
		// TODO(js): Assess the exhaustiveness.
		return nil
	}

	if ndb > 0 {
		switch {
		case ndb == 1 && nhdb == 0:
			a.unsaturation = cmn.UnsaturationDoubleBondC
		case ndb == 1 && nhdb == 1:
			a.unsaturation = cmn.UnsaturationDoubleBondW
		case ndb == 2 && nhdb == 0:
			a.unsaturation = cmn.UnsaturationDoubleBondCC
		case ndb == 2 && nhdb == 1:
			a.unsaturation = cmn.UnsaturationDoubleBondCW
		case ndb == 2 && nhdb == 2:
			a.unsaturation = cmn.UnsaturationDoubleBondWW
		}
		// TODO(js): Assess the exhaustiveness.
	}
	return nil
}

// piElectronCount answers the number of delocalised pi electrons
// contributed by this atom.
//
// This number is important for calculating the aromaticity of the
// rings this atom participates in.
//
// It answers an additional boolean to indicate if the calculation
// could contribute towards computation of aromaticity or not.  A
// `false` value means that the presence of such an atom prevents the
// ring containing it from becoming aromatic.
func (a *_Atom) piElectronCount() (int, bool) {
	mol := a.mol
	wtSum := 100*int16(a.doubleBondCount) + 10*int16(a.singleBondCount) + int16(a.charge)

	switch a.atNum {
	case 6:
		switch wtSum {
		case 19:
			return 2, true
		case 20:
			return 0, true
		case 110:
			return 1, true
		case 120:
			var b *_Bond
			for bid, ok := a.bonds.NextSet(0); ok; bid, ok = a.bonds.NextSet(bid + 1) {
				b = mol.bondWithId(uint16(bid))
				if b.bType == cmn.BondTypeDouble {
					break
				}
			}
			if !b.isCyclic() { // Exocyclic bond.
				return 0, true
			}
			return 1, true // Double bond is in a ring.
		default:
			return 0, true
		}

	case 7:
		switch wtSum {
		case 20, 30:
			return 2, true
		case 110, 121:
			return 1, true
		default:
			return 0, true
		}

	case 8:
		switch wtSum {
		case 20:
			return 2, true
		case 111:
			return 1, true
		default:
			return 0, true
		}

	case 16:
		switch wtSum {
		case 20:
			return 2, true
		case 111:
			return 1, true
		case 120:
			oaid, b := a.firstDoublyBondedNeighbourId()
			oa := mol.atomWithIid(oaid)
			if oa.atNum == 8 && !b.isCyclic() { // Exocyclic bond with an oxygen.
				return 2, true
			}
			return 0, true // Double bond is in a ring.
		case 220:
			c := 0
			for bid, ok := a.bonds.NextSet(0); ok; bid, ok = a.bonds.NextSet(bid + 1) {
				b := mol.bondWithId(uint16(bid))
				if b.bType == cmn.BondTypeDouble {
					if !b.isCyclic() { // Exocyclic bond.
						c++
					}
				}
			}
			if c > 1 {
				return 0, false
			}
			return 0, true
		default:
			return 0, true
		}
	}

	return 0, true
}

// isCyclic answers if this atom participates in at least one ring.
func (a *_Atom) isCyclic() bool {
	return a.rings.Count() > 0
}

// isJunction answers if this atom has more than 2 distinct
// neighbours.
func (a *_Atom) isJunction() bool {
	return a.bonds.Count() > 2
}

// addBond adds the given bond to this atom, if it is not already
// present.  It also adjusts the list of its neighbours appropriately.
//
// Note that it does NOT check to see if the addition conforms to this
// atom's current valence configuration.
func (a *_Atom) addBond(b *_Bond) {
	for bid, ok := a.bonds.NextSet(0); ok; bid, ok = a.bonds.NextSet(bid + 1) {
		if uint16(bid) == b.id {
			return
		}
	}

	a.bonds.Set(uint(b.id))
	nbrId := b.otherAtomIid(a.iId)
	n := int(b.bType)
	for i := 0; i < n; i++ {
		a.nbrs = append(a.nbrs, nbrId)
	}

	switch n {
	case 1:
		a.singleBondCount++
	case 2:
		a.doubleBondCount++
	case 3:
		a.tripleBondCount++
	}
}

// removeBond removes the given bond from this atom, and adjusts the
// list of neighbours appropriately.
//
// Note that it does NOT check to see if the removal conforms to this
// atom's current valence configuration.
func (a *_Atom) removeBond(b *_Bond) {
	nbrId := b.otherAtomIid(a.iId)

	switch b.bType {
	case 1:
		a.singleBondCount--
	case 2:
		a.doubleBondCount--
	case 3:
		a.tripleBondCount--
	}

	wid := 0
	for _, nid := range a.nbrs {
		if nid == nbrId {
			continue
		}
		a.nbrs[wid] = nid
		wid++
	}
	a.nbrs = a.nbrs[:wid]

	a.bonds.Clear(uint(b.id))
}

// bondTo answers the bond that binds this atom to the given atom, if
// one such bond exists.  Answers `nil` otherwise.
func (a *_Atom) bondTo(other uint16) *_Bond {
	mol := a.mol
	for bid, ok := a.bonds.NextSet(0); ok; bid, ok = a.bonds.NextSet(bid + 1) {
		b := mol.bondWithId(uint16(bid))
		if b.otherAtomIid(a.iId) == other {
			return b
		}
	}

	return nil
}

// firstDoublyBondedNeighbourId answers this atom's doubly-bonded
// neighbour having the highest priority.
//
// This method assumes that the molecule is already normalised!
// Calling it on a molecule that has not be normalised yet, leads to
// incorrect results, when more than one double bond exists.
func (a *_Atom) firstDoublyBondedNeighbourId() (uint16, *_Bond) {
	if a.doubleBondCount == 0 {
		return 0, nil
	}

	mol := a.mol
	for bid, ok := a.bonds.NextSet(0); ok; bid, ok = a.bonds.NextSet(bid + 1) {
		b := mol.bondWithId(uint16(bid))
		if b.bType == cmn.BondTypeDouble {
			return b.otherAtomIid(a.iId), b
		}
	}

	panic("Should never be here!")
}

// firstMultiplyBondedNeighbourId answers this atom's doubly-bonded
// neighbour having the highest priority.
//
// This method assumes that the molecule is already normalised!
// Calling it on a molecule that has not be normalised yet, leads to
// incorrect results.
func (a *_Atom) firstMultiplyBondedNeighbourId() (uint16, *_Bond) {
	if a.doubleBondCount == 0 && a.tripleBondCount == 0 {
		return 0, nil
	}

	mol := a.mol
	for bid, ok := a.bonds.NextSet(0); ok; bid, ok = a.bonds.NextSet(bid + 1) {
		b := mol.bondWithId(uint16(bid))
		if b.bType >= cmn.BondTypeDouble {
			return b.otherAtomIid(a.iId), b
		}
	}

	panic("Should never be here!")
}

// isInRingOfSize answers if this atom participates in at least one
// ring of the given size.
func (a *_Atom) isInRingOfSize(n int) bool {
	mol := a.mol
	for rid, ok := a.rings.NextSet(0); ok; rid, ok = a.rings.NextSet(rid + 1) {
		r := mol.ringWithId(uint8(rid))
		if r.size() == n {
			return true
		}
	}

	return false
}

// isInRingLargerThan answers if this atom participates in at least
// one ring that is larger than the given number.
func (a *_Atom) isInRingLargerThan(n int) bool {
	mol := a.mol
	for rid, ok := a.rings.NextSet(0); ok; rid, ok = a.rings.NextSet(rid + 1) {
		r := mol.ringWithId(uint8(rid))
		if r.size() > n {
			return true
		}
	}

	return false
}

// smallestRing answers the smallest unique ring in which this atom
// participates.  If no such unique ring exists, an error is answered.
func (a *_Atom) smallestRing() (uint8, error) {
	if !a.isCyclic() {
		return 0, fmt.Errorf("Atom not cyclic.")
	}

	min := int(math.MaxUint8)
	c := 0
	var ret uint8

	mol := a.mol
	for rid, ok := a.rings.NextSet(0); ok; rid, ok = a.rings.NextSet(rid + 1) {
		r := mol.ringWithId(uint8(rid))
		size := r.size()
		if size == min {
			c++
		} else if size < min {
			ret = uint8(rid)
			c = 1
		}
	}

	if c > 1 {
		return 0, fmt.Errorf("Smallest ring size: %d, number of smallest rings: %d.", min, c)
	}

	return ret, nil
}

// isAromatic answers if this atom is part of an aromatic ring.
//
// Note that the actual aromaticity determination is handled by
// `Ring`.  This method merely answers the set flag.
func (a *_Atom) isAromatic() bool {
	return a.isInAroRing
}

// isInHeteroAromaticRing answers if this atom is part of an aromatic
// ring with at least one hetero atom.
func (a *_Atom) isInHeteroAromaticRing() bool {
	if a.isInAroRing && a.atNum != 6 {
		return true // Simplest case!
	}

	mol := a.mol
	for rid, ok := a.rings.NextSet(0); ok; rid, ok = a.rings.NextSet(rid + 1) {
		r := mol.ringWithId(uint8(rid))
		if r.isHeteroAromatic() {
			return true
		}
	}

	return false
}

// haveCommonRings answers if this atom has at least one ring in
// common with the given atom.
func (a *_Atom) haveCommonRings(aiid uint16) bool {
	other := a.mol.atomWithIid(aiid)

	return a.rings.IntersectionCardinality(other.rings) > 0
}

// isInSameRingsAs answers if this atom participates in exactly the
// same rings as the given atom.
func (a *_Atom) isInSameRingsAs(aiid uint16) bool {
	other := a.mol.atomWithIid(aiid)

	return a.rings.Equal(other.rings)
}

// isInAllRingsOf answers if this atom participates in every ring in
// which the given atom does.
//
// Note that this atom may participate in more rings, as well.
func (a *_Atom) isInAllRingsOf(aiid uint16) bool {
	other := a.mol.atomWithIid(aiid)

	return other.rings.DifferenceCardinality(a.rings) == 0
}

// addRing adds the given ring to the list of this atom's rings.
func (a *_Atom) addRing(r *_Ring) {
	a.rings.Set(uint(r.id))
}

// removeRing removes the given ring from the list of this atom's
// rings.
//
// This method should, usually, not be called.  The composition of a
// ring never changes once it is completed.  However, during a
// reaction, a ring may get broken, leading to its death.  The
// constituent atoms are then notified of that death by calling this
// method.
func (a *_Atom) removeRing(r *_Ring) {
	a.rings.Clear(uint(r.id))
}

// functionalGroup answers the primary feature of this atom, if one is
// present.  Answers `0` otherwise.
func (a *_Atom) functionalGroup() uint16 {
	if len(a.features) == 0 {
		return 0
	}

	return a.features[0]
}

// addFeature adds the given feature to this atom's list of features.
func (a *_Atom) addFeature(fid uint16) {
	a.features = append(a.features, fid)
}

// removeFeature removes the first instance of the given feature from
// this atom's list of features, if it exists in it.
//
// Answers `true` upon a successful removal; `false` otherwise.
func (a *_Atom) removeFeature(fid uint16) bool {
	idx := -1
	for i, f := range a.features {
		if f == fid {
			idx = i
			break
		}
	}

	if idx == -1 {
		return false
	}

	a.features = append(a.features[:idx], a.features[idx+1:]...)
	return true
}

// featureCount answers the number of features substituted on this
// atom.
func (a *_Atom) featureCount() int {
	return len(a.features)
}

// hasFeature answers if this atom has at least one substituent
// matching the given group.
func (a *_Atom) hasFeature(fid uint16) bool {
	for _, f := range a.features {
		if f == fid {
			return true
		}
	}

	return false
}

// isFunctional answers if this atom can play an active role in a
// reaction, in a substituting position.
//
// Note that an atom can yet be a reaction centre without being
// functional.
func (a *_Atom) isFunctional() bool {
	if a.atNum != 6 {
		return true
	}
	if len(a.features) > 0 {
		return true
	}
	if a.unsaturation > cmn.UnsaturationNone {
		return true
	}

	return false
}

// electronWithdrawingNeighbourCount answers the total of unsaturated
// and saturated electron-withdrawing neighbours of this atom.
func (a *_Atom) electronWithdrawingNeighbourCount() int {
	return a.unsatEwNbrCount + a.satEwNbrCount
}

// enolicHydrogenCount answers the number of attached hydrogen atoms,
// if this atom is enolic.
//
// For the attached hydrogen atoms to be enolic, this atom must be a
// carbon, must be saturated, must have at least one
// electron-withdrawing neighbour, and must not be a bridgehead.
func (a *_Atom) enolicHydrogenCount() int {
	if a.atNum != 6 {
		return 0
	}
	if a.unsaturation > cmn.UnsaturationNone {
		return 0
	}
	if a.electronWithdrawingNeighbourCount() == 0 {
		return 0
	}
	if a.isBridgeHead {
		return 0
	}

	return int(a.hCount)
}

// isAtomicLeavingGroup answers if this atom can act as a leaving
// group.
func (a *_Atom) isAtomicLeavingGroup() bool {
	return a.bonds.Count() == 1 && a.atNum != 6
}

// isCH2 answers if this atom is a carbon with exactly two hydrogen
// atoms bound to it.
func (a *_Atom) isCH2() bool {
	return a.atNum == 6 && a.hCount == 2
}

// isCH3 answers if this atom is a carbon with exactly three hydrogen
// atoms bound to it.
func (a *_Atom) isCH3() bool {
	return a.atNum == 6 && a.hCount == 3
}

// isCarbonylC answers if this atom is a carbon with exactly one
// double bond with an oxygen atom.
func (a *_Atom) isCarbonylC() bool {
	if a.atNum != 6 {
		return false
	}
	if a.unsaturation != cmn.UnsaturationDoubleBondW {
		return false
	}

	mol := a.mol
	for bid, ok := a.bonds.NextSet(0); ok; bid, ok = a.bonds.NextSet(bid + 1) {
		b := mol.bondWithId(uint16(bid))
		if b.bType == cmn.BondTypeDouble {
			oaid := b.otherAtomIid(a.iId)
			oa := mol.atomWithIid(oaid)
			if oa.atNum == 8 {
				return true
			}
		}
	}

	return false
}

// isHydroxyl answers if this atom is an oxygen with exactly one
// hydrogen atom bound to it.
func (a *_Atom) isHydroxyl() bool {
	return a.atNum == 8 && a.hCount == 1
}

// isOneOfNOS answers if this atom is any of a nitrogen, oxygen or
// sulfur atoms.
func (a *_Atom) isOneOfNOS() bool {
	switch a.atNum {
	case 7, 8, 16:
		return true
	}

	return false
}

// isOneOfNOPS answers if this atom is any of a nitrogen, oxygen or
// sulfur atoms.
func (a *_Atom) isOneOfNOPS() bool {
	switch a.atNum {
	case 7, 8, 15, 16:
		return true
	}

	return false
}

// isSaturatedC answers if this atom is a carbon with four single
// bonds.
func (a *_Atom) isSaturatedC() bool {
	return a.atNum == 6 && a.unsaturation == cmn.UnsaturationNone
}

// isSaturatedCH2 answers if this atom is a carbon with four single
// bonds, two of which are to hydrogen atoms.
func (a *_Atom) isSaturatedCH2() bool {
	return a.isSaturatedC() && a.hCount == 2
}

// isSaturatedCHavingH answers if this atom is a carbon with four
// single bonds, at least one of which is to a hydrogen atom.
func (a *_Atom) isSaturatedCHavingH() bool {
	return a.isSaturatedC() && a.hCount > 0
}

// isSaturatedHavingH answers if this atom has as many single bonds as
// its current valence configuration, at least one of which is to a
// hydrogen atom.
func (a *_Atom) isSaturatedHavingH() bool {
	return a.unsaturation == cmn.UnsaturationNone && a.hCount > 0
}

// isTerminal answers if this atom has only one neighbour.
func (a *_Atom) isTerminal() bool {
	return a.bonds.Count() == 1
}

// isTerminalHeteroAtom answers if this atom is a hetero atom, having
// only one neighbour.
func (a *_Atom) isTerminalHeteroAtom() bool {
	return a.atNum != 6 && a.bonds.Count() == 1
}

// isTerminalO answers if this atom is an oxygen, having only one
// neighbour.
func (a *_Atom) isTerminalO() bool {
	return a.atNum == 8 && a.bonds.Count() == 1
}

// isTrivalentN answers if this atom is a nitrogen, having current
// valence of 3.
func (a *_Atom) isTrivalentN() bool {
	return a.atNum == 7 && a.valence == 3
}

// isElectronDonating answers if this atom can donate one or more
// electrons.
//
// An atom can donate electrons if it is saturated with its natural
// valence, and has no electron-withdrawing neighbours.
//
// TODO(js): Verify this method's authenticity.
func (a *_Atom) isElectronDonating() bool {
	if a.unsatEwNbrCount > 0 || a.satEwNbrCount > 0 ||
		a.unsaturation != cmn.UnsaturationNone {
		return false
	}

	switch a.atNum {
	case 7, 15:
		return a.bonds.Count() <= 3
	case 8, 16:
		return a.bonds.Count() <= 2
	}

	return false
}

// isHalogen answers if this atom is one of fluorine, chlorine,
// bromine or iodine.
func (a *_Atom) isHalogen() bool {
	switch a.atNum {
	case 9, 17, 35, 53:
		return true
	}

	return false
}

// isNH2orOHorSH answers if this atom is a nitrogen with two attached
// hydrogen atoms, or an oxygen with one attached hydrogen atom, or a
// sulfur with one attached hydrogen atom.
func (a *_Atom) isNH2orOHorSH() bool {
	if a.hCount == 0 || a.unsaturation != cmn.UnsaturationNone {
		return false
	}

	switch a.atNum {
	case 7:
		return a.hCount == 2
	case 8, 16:
		return a.hCount == 1
	}

	return false
}
