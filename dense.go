// Dense solutions.
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

// A DenseSolution is a structure with precalculated points of a solution within
// a given interval. It allows getting the value of any point within that
// interval via interpolation. The zero value is meaningless.
type DenseSolution struct {
	points []Point // This can't be empty. Points are ordered increasingly.
	interp Interpolator
}

// Get the earliest time of the points in the solution.
func (d *DenseSolution) Start() float64 {
	return d.points[0].Time
}

// Get the latest time of the points in the solution.
func (d *DenseSolution) End() float64 {
	return d.points[len(d.points)-1].Time
}

func (d *DenseSolution) search(time float64, min int, max int) []float64 {
	if min == max {
		return d.points[min].Value
	}
	if max-min == 1 {
		switch time {
		case d.points[min].Time:
			return d.points[min].Clone().Value
		case d.points[max].Time:
			return d.points[max].Clone().Value
		}
		return d.interp.FindValue(&d.points[min], &d.points[max], time)
	}
	mid := (max - min) / 2
	midtime := d.points[mid].Time
	switch {
	case time == midtime:
		return d.points[mid].Value
	case time > midtime:
		return d.search(time, mid, max)
	default:
		return d.search(time, min, mid)
	}
}

// Get the value for a point in the solution with a given time between
// `d.Start()` and `d.End()`. If the second return value is false, the point
// was not in the interval and the first return value is meaningless.
func (d *DenseSolution) Get(time float64) ([]float64, bool) {
	if time < d.Start() {
		return d.points[0].Value, false
	}
	if time > d.End() {
		return d.points[len(d.points)-1].Value, false
	}
	return d.search(time, 0, len(d.points)-1), true
}

// Return the list of all explicitly-calculated points. This slice *MUST NOT*
// be modified, nor any of the points it contains.
func (d *DenseSolution) InnerPoints() []Point {
	return d.points
}

// Like CacheSolver.PointCoords but for a DenseSolution.
func (d *DenseSolution) PointCoords() [][]float64 {
	size := len(d.points[0].Value)
	result := make([][]float64, size+1)
	for i := range result {
		result[i] = make([]float64, len(d.points))
	}
	for i, p := range d.points {
		result[0][i] = p.Time
		for j, v := range p.Value {
			result[j+1][i] = v
		}
	}
	return result
}

// Solve an initial value problem in an interval `[start, end]` with the given
// method. The points are calculated eagerly and, when a point is requested not
// given by the method, interpolation is used with the given interpolator.
//
// If the second return value is false, the method didn't get to generate any
// point and the first return value is meaningless. If it's true, the first
// return value contains a solution for an interval that contains at most two
// explicitly-stored points outside the interval requested (the two extremes)
// but which might not contain the whole interval requested if the method
// stopped.
func DenseSolve(
	method Method,
	ivp *IVP,
	start float64,
	end float64,
	interp Interpolator,
) (DenseSolution, bool) {
	forward := method.Forward(ivp)
	backward := method.Backward(ivp)
	points := make([]Point, 0)
	current := &ivp.Start
	pendingAddInitial := true
	var prev Point
	var ok bool

	if end < current.Time { // Special case: end is on the left
		for end < current.Time {
			prev = current.Clone()
			if current, ok = backward.Next(); !ok {
				return DenseSolution{}, false
			}
		}
		if current.Time < end {
			points = append(points, prev)
		}
		points = append(points, current.Clone())
		pendingAddInitial = false
	}
	for start < current.Time { // Move backwards until the start
		if current, ok = backward.Next(); !ok {
			break
		}
		points = append(points, current.Clone())
	}
	size := len(points) // Reverse backwards points
	for i := 0; i < size/2; i++ {
		points[i], points[size-i] = points[size-i], points[i]
	}

	current = &ivp.Start
	if start > current.Time { // Special case: start is on the right
		for start > current.Time {
			prev = current.Clone()
			if current, ok = forward.Next(); !ok {
				return DenseSolution{}, false
			}
		}
		if current.Time > start {
			points = append(points, prev)
		}
		points = append(points, current.Clone())
		pendingAddInitial = false
	}
	if pendingAddInitial {
		points = append(points, ivp.Start)
	}
	for end > current.Time { // Move forwards until the end
		if current, ok = forward.Next(); !ok {
			break
		}
		points = append(points, current.Clone())
	}

	return DenseSolution{points: points, interp: interp}, true
}
