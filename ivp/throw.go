// Parabolic throw problem with air friction.
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

// Gravity acceleration in the Earth surface
const EARTH_GRAVITY float64 = 9.806

// A ParabolicThrow is the result of throwing an object, which we consider to
// be punctual, from a given height and with some initial velocity. We assume
// that Newton mechanics apply and that air resistance has a force proportional
// to the square of the velocity, and opposite to it.
//
// The resulting ODE is `x'' = -gj - (k/m)*x'*|x'|`, where `x` is the position,
// `g` is the gravity acceleration, `j` is the vertical unit vector (upwards),
// `k` is the air resistance, and `m` is the mass of the object.
//
// The comments on fields assume standard SI units, but other units may be used.
type ParabolicThrow struct {
	Height     float64    // Initial height over the floor (m).
	Mass       float64    // Mass of the object (kg). *Must* be positive.
	V0         [2]float64 // Initial velocity (x,y) (m/s).
	Resistance float64    // Air resistance (kg/m).
	Gravity    float64    // Gravity acceleration (m/s^2).
}

// Create an IVP from the data on this problem.
//
// Result values have the form (x,y,vx,vy), where (x,y) is the position and
// (vx,vy) is the velocity. The initial x and time are 0.
func (p *ParabolicThrow) ToIVP() mned.IVP {
	ratio := p.Resistance / p.Mass
	gravity := p.Gravity
	return mned.IVP{
		Derivative: func(p mned.Point) ([]float64, bool) {
			x := p.Value
			air := -ratio * math.Sqrt(x[2]*x[2]+x[3]*x[3])
			return []float64{
				x[2], x[3], air * x[2], air*x[3] - gravity,
			}, true
		},
		Start: mned.Point{
			Time:  0,
			Value: []float64{0, p.Height, p.V0[0], p.V0[1]},
		},
	}
}
