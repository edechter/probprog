# definitely #
* undoable effects
* save the number of accepts and rejects

# maybe #
* make pramb a macro like amb so it doesn't evaluate its arguments?

# next version #
* some lazy-evaluating magic for gaussians, conjugacy, etc
    - force should be generic. but then maybe we don't need sample:val (since we
    can extract it using force). yeah the interface should be reworked.
* data is an environment
* special e.g. mean function (which can be exact or smart) with a default
  implementation. generics!
* stream of samples (in terms of resume, should be easy, then function to
  estimate means and stuff should be cleaner)
* everything should be uncollapsed by default and not add to choice chain

# notes #
* generic dispatch doesn't work like it says it does, so uncollapsed gaussian
  samples need to come first in the argument list to +
