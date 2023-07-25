// Example comparing the formula and numeric approximation for Hooke's law
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
const tol = 0.01

var methods = [4]methodDescriptor{
	{name: "Euler", m: method.Euler(step)},
	{name: "RK4", m: method.RK4(step)},
	{name: "Adaptive RK4", m: method.AdaptiveRK4(tol, step, 1e-6, 1)},
	{name: "RK-Fehlberg", m: method.RKFehlberg(tol, step, 1e-6, 1)},
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
	spring := ivp.HookeSpring{
		Mass:   1,
		Spring: 0.7,
		X0:     1,
		Length: 0.7,
	}
	problem := spring.ToIVP()
	amplitude := spring.X0 - spring.Length
	sqratio := math.Sqrt(spring.Spring / spring.Mass)
	realSol := func(t float64) []float64 {
		return []float64{
			amplitude*math.Cos(sqratio*t) + spring.Length,
			-amplitude * sqratio * math.Sin(sqratio*t),
		}
	}
	interp := mned.HermiteForIVP(&problem)

	for _, m := range methods {
		solution, _ := mned.DenseSolve(
			m.m, &problem, 0, 40, interp,
		)
		points := solution.InnerPoints()
		fmt.Printf("%v: Error of %v in %v steps.\n",
			m.name, maxerr(points, realSol), len(points))
	}
}
