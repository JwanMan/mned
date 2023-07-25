// Spring movement problem based on Hooke's law.
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

package ivp

import (
	"math"
	"github.com/JwanMan/mned"
)

// A HookeSpring is a spring moving following Hooke's law.
//
// The spring is laid out horizontally over a table and has one end attached to
// a body which receives a certain friction from the table and such that the
// mass of the spring is negligible compared to the mass of the body. The other
// end of the spring is attached to a harmonic oscillator. The position of the
// body is measured as the (signed) distance to the rest position of the
// oscillator, and the velocity is measured in terms of such position.
//
// Only Mass and Spring are required, although the result of only specifying
// that would be really boring. Note also that, since we're following Hooke's
// law, the length of the spring can get negative.
type HookeSpring struct {
	Mass       float64 // Mass of the body (kg).
	Spring     float64 // Hooke's constant of the spring (kg/s^2).
	X0         float64 // Initial position of the body (m).
	V0         float64 // Initial velocity of the body (m/s).
	Length     float64 // Length of the spring at rest (m).
	Friction   float64 // Friction between the body and the table (kg/s).
	Amplitude  float64 // Amplitude of the oscillator's force (N).
	AngleSpeed float64 // Angle speed of the oscillator (rad/s).
}

// Get the IVP associated to the problem. The solution values have the form
// (length, speed), time starts at 0, and the harmonic oscillator starts at
// its natural position.
//
// The ODE is `m * x''(t) = A*sin(wt) - k*(x(t)-l) - rx'(t)`, where m is the
// mass of the body, x is its position, A is the amplitude of the oscillator,
// w is its angular velocity, k is the spring's constant, l is its length at
// rest, and r is the friction coefficient.
func (h *HookeSpring) ToIVP() mned.IVP {
	adjSpring := h.Spring / h.Mass
	adjFriction := h.Friction / h.Mass
	adjAmplitude := h.Amplitude / h.Mass
	w := h.AngleSpeed
	length := h.Length

	return mned.IVP{
		Derivative: func(p mned.Point) ([]float64, bool) {
			return []float64{
				p.Value[1],
				adjAmplitude*math.Sin(w*p.Time) -
					adjSpring*(p.Value[0]-length) -
					adjFriction*p.Value[1],
			}, true
		},
		Start: mned.Point{
			Time:  0,
			Value: []float64{h.X0, h.V0},
		},
	}
}
