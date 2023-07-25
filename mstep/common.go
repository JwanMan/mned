// Common code to define stepping methods.
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

type Step struct {
	Value []float64
	Deriv []float64
}

func Implicit(steps []Step, a []float64, b []float64, step float64) []float64 {
	if len(steps) == 0 {
		panic("There must be at least one step!\n")
	}
	if len(steps) < len(a) || len(steps) < len(b) {
		panic("There must be at least as many steps as coefficients.\n")
	}

	result := make([]float64, len(steps[0].Value))
	for i, v := range a {
		for j := range result {
			result[j] += v * steps[i].Value[j]
		}
	}
	for i, v := range b {
		for j := range result {
			result[j] += step * v * steps[i].Deriv[j]
		}
	}
	return result
}

func Shift(steps []Step) {
	for i := len(steps) - 1; i > 0; i-- {
		steps[i] = steps[i-1]
	}
}
