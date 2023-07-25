// Example visualization of Hooke's law.
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
	"mned"
	"mned/ivp"
	"mned/method"
)

func spring1(speed float64) {
	spring := ivp.HookeSpring{
		Mass:       1,
		Spring:     1.5,
		Length:     0.7,
		Friction:   0.3,
		Amplitude:  0.4,
		AngleSpeed: speed,
	}
	problem := spring.ToIVP()
	plot, err := glot.NewPlot(2, true, false)
	if err != nil {
		fmt.Printf("Creating the plot: %v\n", err)
		return
	}
	if err := plot.SetXLabel(
		fmt.Sprintf("Time (s) -- w = %v", speed),
	); err != nil {
		fmt.Printf("Drawing the X label: %v\n", err)
		return
	}

	solution, _ := mned.DenseSolve(
		method.RKFehlberg(0.0001, 0.01, 1e-7, 0.5),
		&problem,
		0,
		40,
		mned.HermiteForIVP(&problem),
	)
	coord := solution.PointCoords()
	if err := plot.AddPointGroup("Position (m)", "lines", [][]float64{
		coord[0], coord[1],
	}); err != nil {
		fmt.Printf("Drawing the position: %v\n", err)
	}
	if err := plot.AddPointGroup("Velocity (m/s)", "lines", [][]float64{
		coord[0], coord[2],
	}); err != nil {
		fmt.Printf("Drawing the velocity: %v\n", err)
	}
}

func main() {
	spring1(0)
	spring1(1.3)
	spring1(2.4)
	spring1(3.5)
}
