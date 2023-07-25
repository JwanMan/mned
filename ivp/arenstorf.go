// Arenstorf orbit problem.
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

var arenstorfOrbit = mned.IVP{
	Derivative: func(p mned.Point) ([]float64, bool) {
		const μ = 0.012277471
		const μ2 = 1 - μ
		xμ := p.Value[0] + μ
		xμ2 := p.Value[0] - μ2
		D1 := xμ*xμ + p.Value[1]*p.Value[1]
		D1 *= math.Sqrt(D1)
		D2 := xμ2*xμ2 + p.Value[1]*p.Value[1]
		D2 *= math.Sqrt(D2)
		return []float64{
			p.Value[2],
			p.Value[3],
			p.Value[0] + 2*p.Value[3] - μ2*xμ/D1 - μ*xμ2/D2,
			p.Value[1]*(1-μ2/D1-μ/D2) - 2*p.Value[2],
		}, true
	},
	Start: mned.Point{
		Time:  0,
		Value: []float64{0.994, 0, 0, -2.001585106},
	},
}

// The Arenstorf orbit is a stable orbit of a satellite around the Earth taking
// into account the gravity of the Earth and the Moon and considering the
// satellite to have a negligible mass compared to both. The orbit was
// discovered by Richard Arenstorf and is given by the following formulae:
//
// ```
// μ = 0.012277471, μ' = 1 - μ;
// D1 = ((x+μ)^2 + y^2)^(3/2), D2 = ((x-μ')^2 + y^2)^(3/2);
// x" = x + 2y' - μ'(x+μ)/D1 - μ(x-μ')/D2;
// y" = y - 2x' - μ'y/D1 - μy/D2;
// x(0) = 0.994, y(0) = 0, x'(0) = 0, y'(0) = -2.001585106.
// ```
//
// The orbit has three big loops at one side and a small one at the other, and
// varying μ changes the number of loops. Because of the unstability, this
// problem is often used to test different ODE solving methods.
func ArenstorfOrbit() mned.IVP {
	return arenstorfOrbit.Clone()
}
