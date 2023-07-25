// Events mechanism for early iteration ending.
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

package mned

import "math"

// An Event is a condition that can happen in a point in the solution of an
// initial value problem and an associated action to take. The event happens
// at the point of the solution where a given function changes its sign.
//
// The Cross function is a continuous function whose zeroes are the points
// where the event happens. The Tolerance is the margin of error allowed; the
// maximum absolute value of `Cross(p)` such that `p` is considered to be close
// enough to an event to be passed to Action. The Action is to be called when
// an occurrence of the event is found; it can use the point in some way and it
// returns a boolean indicating whether the calculation should continue.
type Event struct {
	Cross     func(*Point) float64
	Tolerance float64
	Action    func(*Point) bool
}

// If e.Cross has a zero between p1 and p2, find such a zero or a point `p`
// close enough to the zero that `|e.Cross(p)| < e.Tolerance`. To get the
// intermediante points, the given interpolator is used.
func (e *Event) FindPoint(i Interpolator, p1 *Point, p2 *Point) Point {
	var min, max, mid Point
	if p1.Time < p2.Time {
		min, max = *p1, *p2
	} else {
		max, min = *p1, *p2
	}
	downwards := e.Cross(p1) > 0
	for {
		mid.Time = (min.Time + max.Time) / 2
		mid.Value = i.FindValue(p1, p2, mid.Time)
		value := e.Cross(&mid)
		if math.Abs(value) < e.Tolerance {
			return mid
		}
		if downwards {
			if value > 0 {
				min = mid
			} else {
				max = mid
			}
		} else {
			if value > 0 {
				max = mid
			} else {
				min = mid
			}
		}
	}
}

type requiredAction struct {
	index int
	point Point
}

type requiredActions []requiredAction

func (r requiredActions) Len() int {
	return len(r)
}

func (r requiredActions) Less(i, j int) bool {
	return r[i].point.Time < r[j].point.Time
}

func (r requiredActions) Swap(i, j int) {
	r[i], r[j] = r[j], r[i]
}

func differentSign(f1 float64, f2 float64) bool {
	return (f1 <= 0 && f2 >= 0) || (f1 >= 0 && f2 <= 0)
}
