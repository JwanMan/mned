// Example visualization of parabolic throw with air friction.
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

func addPlot(
	plot *glot.Plot, name string, arrays [][]float64, idx int, idy int,
) {
	if err := plot.AddPointGroup(
		name,
		"lines",
		[][]float64{arrays[idx], arrays[idy]}); err != nil {
		fmt.Printf("Error adding %v: %w\n", name, err)
		os.Exit(1)
	}
}

func main() {
	var end mned.Point
	var err error

	throw := ivp.ParabolicThrow{
		Height:     300,
		V0:         [2]float64{100, 0},
		Mass:       1,
		Resistance: 0.01,
		Gravity:    ivp.EARTH_GRAVITY,
	}
	problem := throw.ToIVP()
	solution := mned.CacheSolve(
		method.Euler(0.01),
		&problem,
		mned.LinearInterpolator{},
		mned.Event{
			Cross: func(p *mned.Point) float64 {
				return p.Value[1]
			},
			Tolerance: 0.01,
			Action: func(p *mned.Point) bool {
				end = *p
				return false
			},
		})
	solution.StepToEnd()
	fmt.Printf(
		"Reached ground after %v s at %v m, with (%v, %v) m/s.\n",
		end.Time,
		end.Value[0],
		end.Value[2],
		end.Value[3],
	)

	arrays := solution.PointCoords()
	byTime, err := glot.NewPlot(2, true, false)
	if err != nil {
		fmt.Printf("Creating plot window: %w\n")
		os.Exit(1)
	}
	defer byTime.Close()
	if err := byTime.SetXLabel("Time (s)"); err != nil {
		fmt.Printf("Adding X label: %w\n")
		os.Exit(1)
	}
	addPlot(byTime, "X position (m)", arrays, 0, 1)
	addPlot(byTime, "Y position (m)", arrays, 0, 2)
	addPlot(byTime, "X speed (m/s)", arrays, 0, 3)
	addPlot(byTime, "Y speed (m/s)", arrays, 0, 4)

	byTraj, err := glot.NewPlot(2, true, false)
	if err != nil {
		fmt.Printf("Creating plot window: %w\n")
		os.Exit(1)
	}
	defer byTraj.Close()
	if err := byTraj.SetXLabel("X coordinate"); err != nil {
		fmt.Printf("Adding X label: %w\n")
		os.Exit(1)
	}
	if err := byTraj.SetYLabel("Y coordinate"); err != nil {
		fmt.Printf("Adding Y label: %w\n")
		os.Exit(1)
	}
	addPlot(byTraj, "Trajectory (m)", arrays, 1, 2)
	addPlot(byTraj, "Speed trajectory (m/s)", arrays, 3, 4)
}
