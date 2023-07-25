// Interpolation routines.
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

// An Interpolator tries to approximate a function based on a finite set of
// known points in it. We'll focus just on interpolation with two points since
// this is precise enough for our purposes, although the interpolator might be
// constructed with some knowledge about the underlying function.
type Interpolator interface {
	// Approximate the value of a function in a point `t` given two
	// **different** points `p1` and `p2` of the function such that `t`
	// is between `p1.Time` and `p2.Time`.
	FindValue(p1 *Point, p2 *Point, t float64) []float64
}

// A LinearInterpolator is an interpolator that assumes that there's a straight
// line between the two points given. It can be constructed with `new`.
type LinearInterpolator struct{}

func (li LinearInterpolator) FindValue(
	p1 *Point, p2 *Point, t float64,
) []float64 {
	// The line between (t1,x1) and (t2,x2) is x1 + (x2-x1)*(t-t1)/(t2-t1).
	// If ratio:=(t-t1)/(t2-t1), this is x2*ratio + x1*(1-ratio)
	size := len(p1.Value)
	result := make([]float64, size)
	ratio := (t - p1.Time) / (p2.Time - p1.Time)
	for i := 0; i < size; i++ {
		result[i] = p1.Value[i]*(1-ratio) + p2.Value[i]*ratio
	}
	return result
}

// A HermiteInterpolator is an interpolator that takes into account the
// derivative of the solution function in the end points.
type HermiteInterpolator struct {
	F func(Point) ([]float64, bool) // The derivative function.
}

func (h HermiteInterpolator) FindValue(
	p1 *Point, p2 *Point, t float64,
) []float64 {
	diff := p2.Time - p1.Time
	off := t - p1.Time
	d01, ok := h.F(*p1)
	if !ok { // Can't derive: fallback
		return LinearInterpolator{}.FindValue(p1, p2, t)
	}
	d23, ok := h.F(*p2)
	if !ok {
		return LinearInterpolator{}.FindValue(p1, p2, t)
	}
	d12 := make([]float64, len(d01))
	for i := range d12 {
		d12[i] = (p2.Value[i] - p1.Value[i]) / diff
	}
	d02 := make([]float64, len(d01))
	for i := range d02 {
		d02[i] = (d12[i] - d01[i]) / diff
	}
	d13 := make([]float64, len(d01))
	for i := range d13 {
		d13[i] = (d23[i] - d12[i]) / diff
	}
	d03 := make([]float64, len(d01))
	for i := range d03 {
		d03[i] = (d13[i] - d02[i]) / diff
	}
	result := make([]float64, len(d01))
	for i := range result {
		result[i] = p1.Value[i] + off*
			(d01[i]+off*(d02[i]+(off-diff)*d03[i]))
	}
	return result
}

// Create a Hermite interpolator suitable for the given IVP.
func HermiteForIVP(ivp *IVP) HermiteInterpolator {
	return HermiteInterpolator{F: ivp.Derivative}
}
