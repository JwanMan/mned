// Example comparing the formula and numeric approximation for parabolic throw.
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

package main

import (
	"fmt"
	"math"
	"github.com/JwanMan/mned"
	"github.com/JwanMan/mned/ivp"
	"github.com/JwanMan/mned/method"
)

type methodDescriptor struct {
	name string
	m    mned.Method
}

const step = 0.01

var methods = [3]methodDescriptor{
	{name: "Euler", m: method.Euler(step)},
	{name: "Modified Euler", m: method.ModifiedEuler(step)},
	{name: "RK4", m: method.RK4(step)},
}

func normdiff(p1 []float64, p2 []float64) float64 {
	var norm float64 = 0
	for i, v := range p1 {
		norm += (v - p2[i]) * (v - p2[i])
	}
	return norm
}

func maxerr(points []mned.Point, sol func(float64) []float64) float64 {
	var max float64 = 0
	for _, p := range points {
		actual := sol(p.Time)
		norm := normdiff(actual, p.Value)
		if norm > max {
			max = norm
		}
	}
	return max
}

func main() {
	// Here we test the precision of different methods against the known
	// problem of parabolic throw without air resistance. The problem is
	// x''(t) = -gj, so x(t) = -g(t^2/2)j + v0*t + hj, thus
	// x(t)[1] = 0 when t = (v0[1] + sqrt(v0[1]^2 + 2hg)) / g
	var fallen float64
	throw := ivp.ParabolicThrow{
		Height:  300,
		V0:      [2]float64{100, 0},
		Mass:    1,
		Gravity: ivp.EARTH_GRAVITY,
	}
	problem := throw.ToIVP()
	realSolution := func(t float64) []float64 {
		return []float64{
			throw.V0[0] * t,
			throw.Height + throw.V0[1]*t - throw.Gravity*t*t/2,
			throw.V0[0],
			throw.V0[1] - throw.Gravity*t,
		}
	}
	realEnd := (throw.V0[1] + math.Sqrt(
		throw.V0[1]*throw.V0[1]+2*throw.Height*throw.Gravity,
	)) / throw.Gravity * throw.V0[0]

	for _, method := range methods {
		solution := mned.CacheSolve(
			method.m, &problem, mned.LinearInterpolator{},
			mned.Event{
				Cross: func(p *mned.Point) float64 {
					return p.Value[1]
				},
				Tolerance: 0.01,
				Action: func(p *mned.Point) bool {
					fallen = p.Value[0]
					return false
				},
			})

		solution.StepToEnd()
		soldiff := maxerr(solution.ForwardPoints(), realSolution)
		enddiff := math.Abs(fallen - realEnd)
		fmt.Printf("%v: Point diff @ %v, end diff @ %v.\n",
			method.name, soldiff, enddiff)
	}
}
