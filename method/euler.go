// Euler method and variants.
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

func eulerStep(ivp *mned.IVP, h float64) bool {
	deriv, ok := ivp.Derivative(ivp.Start)
	if !ok {
		return false
	}
	for i := 0; i < len(deriv); i++ {
		ivp.Start.Value[i] += h * deriv[i]
	}
	return true
}

// Make an instance of the Euler method with the given step, which must be
// positive. For a problem `x'(t)=f(t,x(t)), x(t0)=x0`, a step of the Euler
// method is given by `x(t0+step) = x(t0) + step*f(t0,x(t0))`.
func Euler(step float64) mned.Method {
	return mned.FixedStepMethod{Step: step, Method: eulerStep}
}

// Make an instance of the adaptive version of the Euler method given the
// tolerance, the initial step, and the minimum and maximum steps.
func AdaptiveEuler(
	tolerance float64, step float64, hmin float64, hmax float64,
) mned.Method {
	return &mned.AdaptiveStepMethod{
		Step:      step,
		Hmin:      hmin,
		Hmax:      hmax,
		Tolerance: tolerance,
		Method:    eulerStep,
		Order:     1,
	}
}

func modEulerStep(ivp *mned.IVP, h float64) bool {
	deriv, ok := ivp.Derivative(ivp.Start)
	if !ok {
		return false
	}
	point := ivp.Start.Clone()
	point.Time += h
	for i, v := range deriv {
		point.Value[i] += h * v
	}
	deriv2, ok := ivp.Derivative(point)
	for i, _ := range deriv {
		ivp.Start.Value[i] += h * (deriv[i] + deriv2[i]) / 2
	}
	return true
}

// Make an instance of the modified Euler method with the given step, which must
// be positive. This is a FixedStepMethod given by
// `x(t+h)=x(t) + h/2 * (f(t,x(t)) + f(t+h,x(t)+h*f(t,x(t))))`.
func ModifiedEuler(step float64) mned.Method {
	return mned.FixedStepMethod{Step: step, Method: modEulerStep}
}
