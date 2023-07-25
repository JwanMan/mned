// Example of trapezium method for exponential solutions.
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
	"github.com/JwanMan/mned/method"
	"github.com/JwanMan/mned/mstep"
)

type methodDescriptor struct {
	name string
	m    mned.Method
}

const step = 0.01
const tol = 0.01

var df mstep.YDerivative = func(t float64, v float64) float64 { return -1 }
var methods = [5]methodDescriptor{
	{name: "Euler", m: method.Euler(0.01)},
	{name: "RK4", m: method.RK4(0.1)},
	{name: "Implicit Euler", m: mstep.ImplicitEuler(0.01, 0.01, df)},
	{name: "Trapezium", m: mstep.Trapezium(0.01, 0.01, df)},
	{name: "RKF", m: method.RKFehlberg(0.01, 0.01, 1e-7, 0.1)},
}

func main() {
	problem := mned.IVP{
		Derivative: func(p mned.Point) ([]float64, bool) {
			return []float64{-p.Value[0]}, true
		},
		Start: mned.Point{Time: 0, Value: []float64{1}},
	}
	for _, m := range methods {
		solution, _ := mned.DenseSolve(
			m.m, &problem, 0, 10, mned.HermiteForIVP(&problem),
		)
		points := solution.InnerPoints()
		maxerr := float64(0)
		for _, p := range points {
			realVal := math.Exp(-p.Time)
			err := math.Abs(p.Value[0] / realVal)
			if err < 1 {
				err = 1 / err
			}
			err -= 1
			if err > maxerr {
				maxerr = err
			}
		}
		fmt.Printf("%v method got %v steps, max rel err = %v.\n",
			m.name, len(points), maxerr)
	}
}
