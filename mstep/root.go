// Root-finding algorithms.
//
// Copyright (C) 2020 Juan Marín Noguera
//
// This file is part of Solvned.
//
// Solvned is free software: you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option) any 
// later version.
//
// Solvned is distributed in the hope that it will be useful, but WITHOUT ANY 
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Solvned. If not, see <https://www.gnu.org/licenses/>.

package mstep

import "math"

func NewtonRoot1D(
	x0, tolerance float64, f, df func(float64) (float64, bool),
) (float64, bool) {
	prev := x0
	for {
		fprev, ok := f(prev)
		if !ok {
			return 0, false
		}
		dfprev, ok := df(prev)
		if !ok {
			return 0, false
		}
		off := fprev / dfprev
		next := prev - off
		if math.Abs(off) <= tolerance {
			return next, true
		}
		prev = next
	}
}

func SecantRoot1D(
	x0, x1, tolerance float64, f func(float64) (float64, bool),
) (float64, bool) {
	prev := x0
	fprev, ok := f(prev)
	if !ok {
		return 0, false
	}
	cur := x1
	fcur, ok := f(cur)
	if !ok {
		return 0, false
	}
	for {
		off := fcur * (prev - cur) / (fprev - fcur)
		next := cur - off
		if math.Abs(off) <= tolerance {
			return next, true
		}
		fnext, ok := f(next)
		if !ok {
			return 0, false
		}
		prev, fprev, cur, fcur = cur, fcur, next, fnext
	}
}
