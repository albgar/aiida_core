Suggestions for Siesta flow rethink

The plugin work has uncovered a few issues worth considering

* Essential and non-essential output files

Just one "structure" file should be needed...

* Appropriate places for output

Forces and Stresses should be output (in CML also) at each "state_analysis"
step. Spin should too (it is currently in "siesta_analysis" at the
end).

* Restarting capabilities

What should be the specifications?

* Basis set specification

It cannot be made too simple, but it can be streamlined.

