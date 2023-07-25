// Runge-Kutta method and adaptive variant.
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

package method

import "github.com/JwanMan/mned"

func RK4Step(ivp *mned.IVP, h float64) ([]float64, bool) {
	k1, ok := ivp.Derivative(ivp.Start)
	if !ok {
		return nil, false
	}
	for i, _ := range k1 {
		k1[i] *= h
	}

	p1 := ivp.Start.Clone()
	p1.Time += h / 2
	for i, v := range k1 {
		p1.Value[i] += v / 2
	}
	k2, ok := ivp.Derivative(p1)
	if !ok {
		return nil, false
	}
	for i, _ := range k2 {
		k2[i] *= h
	}

	p2 := ivp.Start.Clone()
	p2.Time += h / 2
	for i, v := range k2 {
		p2.Value[i] += v / 2
	}
	k3, ok := ivp.Derivative(p2)
	if !ok {
		return nil, false
	}
	for i, _ := range k3 {
		k3[i] *= h
	}

	p3 := ivp.Start.Clone()
	p3.Time += h
	for i, v := range k3 {
		p3.Value[i] += v
	}
	k4, ok := ivp.Derivative(p3)
	if !ok {
		return nil, false
	}
	for i, _ := range k4 {
		k4[i] *= h
	}

	for i, _ := range ivp.Start.Value {
		ivp.Start.Value[i] += (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) / 6
	}
	return k1, true
}

func rk4Step(ivp *mned.IVP, h float64) bool {
	_, ok := RK4Step(ivp, h)
	return ok
}

// Initialize the 4th-order Runge-Kutta method, a.k.a. RK4. This is a fixed step
// method given by the following formulae:
// ```
// k1 = h*f(t, x(t))
// k2 = h*f(t+h/2, x(t)+k1/2)
// k3 = h*f(t+h/2, x(t)+k2/2)
// k4 = h*f(t+h, x(t)+k3)
// x(t+h) = x(t) + (k1 + 2*k2 + 2*k3 + k4) / 6
// ```
func RK4(step float64) mned.Method {
	return mned.FixedStepMethod{Step: step, Method: rk4Step}
}

func AdaptiveRK4(
	tolerance float64, step float64, hmin float64, hmax float64,
) mned.Method {
	return &mned.AdaptiveStepMethod{
		Step:      step,
		Hmin:      hmin,
		Hmax:      hmax,
		Tolerance: tolerance,
		Method:    rk4Step,
		Order:     4,
	}
}
