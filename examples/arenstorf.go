// Example visualization of the Arenstorf orbit.
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
	"github.com/Arafatk/glot"
	"github.com/JwanMan/mned"
	"github.com/JwanMan/mned/ivp"
	"github.com/JwanMan/mned/method"
	"os"
)

func main() {
	var backed bool = false
	problem := ivp.ArenstorfOrbit()
	solution := mned.CacheSolve(
		method.RKFehlberg(0.01, 0.01, 1e-7, 0.1),
		&problem,
		mned.HermiteForIVP(&problem),
		mned.Event{
			Cross: func(p *mned.Point) float64 {
				return p.Value[1]
			},
			Tolerance: 0.01,
			Action: func(p *mned.Point) bool {
				if p.Value[0] < 0 {
					backed = true
				}
				return true
			},
		},
		mned.Event{
			Cross: func(p *mned.Point) float64 {
				return p.Value[1]
			},
			Tolerance: 0.01,
			Action: func(p *mned.Point) bool {
				return p.Value[0] <= 0.9 || !backed
			},
		},
	)
	solution.StepToEnd()
	coords := solution.PointCoords()
	fmt.Printf("Took %v steps.\n", len(coords[0]))
	plot, err := glot.NewPlot(2, true, false)
	if err != nil {
		fmt.Printf("Couldn't create plot: %v\n", err)
		os.Exit(1)
	}
	if err := plot.AddPointGroup("Orbit", "lines", [][]float64{
		coords[1], coords[2],
	}); err != nil {
		fmt.Printf("Couldn't draw the points: %v\n", err)
		os.Exit(1)
	}
}
