# besmarts-smee-refit

This is a WIP repo for working on an experimental force field split and re-fit using BeSMARTS and Smee.

## Problem description

After going through Sage 2.2.1's performance, several torsions have been identified with our existing torsions that could potentially be fixed with some combination of splitting torsions, modifying the torsional form, and re-fitting.

* t48a: this torsion is only fit to one molecule where the complementary torsion drives most of the energy. The term may benefit from splitting to cover aromatic and non-aromatic systems, as well as a change in the torsional term
* t17: Across a symmetric ring where the carbon adjacent to the ring has one additional substituent (e.g. ethyl benzene) the Sage energy is very high (mostly contributed by nonbonded and angular terms). The torsional potential should possibly have a different functional form that mitigates the high MM energy. The term is also maybe too generic, as the profiles of higher order substitution yield a better match.
* t19, which is mostly trained to terminal methyls, has an odd functional form where the periodicity=1 term dominates and gives each single hydrogen a strong central peak instead of weak n=3 peaks
* t18 performs poorly for amide-adjacent torsions
* t105 covers an O linker with an sp2 or sp3 terminus. While the sp3 profiles match the QM well, the sp2 profiles look too stiff

## Objective

Using smee as the force field fitting and evaluation backend and BeSmarts to generate potential clusters or splits ("split/merge") and torsion functional forms ("modify"). Data would comprise optimization *and* torsiondrive data.

## Data

Data is curated from both existing training data (specifically the Sage 2.2.1 set) as well as the additional molecules suggested by Cresset.


## Some notes on parameters

Cresset included some ideas as to what may be the correct answer based on their analysis.

* t48a: remove the periodicity=3 term and possibly split between aromatic and non-aromatic systems.
    * potentially: split out t47 with an O-specific type for HCCO, as currently t47 may be too strong
* t17: either split out the term or have an additional 2-fold term
* t19: a good end state would result in a stronger n=3 term and maybe even removal of the n=1
* t105 may need splitting