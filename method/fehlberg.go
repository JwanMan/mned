// Runge-Kutta-Fehlberg method.
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

import (
	"math"
	"github.com/JwanMan/mned"
)

type rkFehlbergStepper struct {
	current mned.Point
	deriv   func(mned.Point) ([]float64, bool)
	step    float64
	tol     float64
	hmin    float64
	hmax    float64
}

func (s *rkFehlbergStepper) Next() (*mned.Point, bool) {
	for s.step >= s.hmin {
		k1, ok := s.deriv(s.current)
		if !ok {
			s.step *= 0.5
			continue
		}
		for i := range k1 {
			k1[i] *= s.step
		}

		p2 := s.current.Clone()
		p2.Time += s.step / 4
		for i := range p2.Value {
			p2.Value[i] += k1[i] / 4
		}
		k2, ok := s.deriv(p2)
		if !ok {
			s.step *= 0.5
			continue
		}
		for i := range k2 {
			k2[i] *= s.step
		}

		p3 := s.current.Clone()
		p3.Time += s.step * 3 / 8
		for i := range p3.Value {
			p3.Value[i] += (3*k1[i] + 9*k2[i]) / 32
		}
		k3, ok := s.deriv(p3)
		if !ok {
			s.step *= 0.5
			continue
		}
		for i := range k3 {
			k3[i] *= s.step
		}

		p4 := s.current.Clone()
		p4.Time += s.step * 12 / 13
		for i := range p4.Value {
			p4.Value[i] +=
				(1932*k1[i] - 7200*k2[i] + 7296*k3[i]) / 2197
		}
		k4, ok := s.deriv(p4)
		if !ok {
			s.step *= 0.5
			continue
		}
		for i := range k4 {
			k4[i] *= s.step
		}

		p5 := s.current.Clone()
		p5.Time += s.step
		for i := range p5.Value {
			p5.Value[i] += 439*k1[i]/216 - 8*k2[i] +
				3680*k3[i]/513 - 845*k4[i]/4104
		}
		k5, ok := s.deriv(p5)
		if !ok {
			s.step *= 0.5
			continue
		}
		for i := range k5 {
			k5[i] *= s.step
		}

		p6 := s.current.Clone()
		p6.Time += s.step / 2
		for i := range p6.Value {
			p6.Value[i] +=
				-8*k1[i]/27 + 2*k2[i] - 3544*k3[i]/2656 +
					1859*k4[i]/4104 - 11*k5[i]/40
		}
		k6, ok := s.deriv(p6)
		if !ok {
			s.step *= 0.5
			continue
		}
		for i := range k6 {
			k6[i] *= s.step
		}

		sqerr := float64(0)
		// Get \tilde{w}_{i+1}-w_{i+1} directly:
		// * (- 16/135 25/216)
		// 1/360
		// * (- 6656/12825 1408/2565)
		// -128/4275
		// * (- 28561/56430 2197/4104)
		// -2197/75240
		// * (- -9/50 -1/5)
		// 1/50
		for i := range s.current.Value {
			v := k1[i]/360 - 128*k3[i]/4275 -
				2197*k4[i]/75240 + k5[i]/50 + 2*k6[i]/55
			sqerr += v * v
		}
		err := math.Sqrt(sqerr)

		if err > s.tol {
			s.step = mned.AdaptiveUpdateStep(
				s.step, s.tol, err, s.hmax, 4,
			)
			continue
		}

		for i := range s.current.Value {
			s.current.Value[i] += 25*k1[i]/216 + 1408*k3[i]/2565 +
				2197*k4[i]/4104 - k5[i]/5
		}
		s.current.Time += s.step
		s.step = mned.AdaptiveUpdateStep(s.step, s.tol, err, s.hmax, 4)
		return &s.current, true
	}
	return nil, false
}

type rkFehlbergMethod struct {
	step float64
	tol  float64
	hmin float64
	hmax float64
}

func (m *rkFehlbergMethod) withStep(step float64, ivp *mned.IVP) mned.Stepper {
	return &rkFehlbergStepper{
		step:    step,
		tol:     m.tol,
		hmin:    m.hmin,
		hmax:    m.hmax,
		current: ivp.Start.Clone(),
		deriv:   ivp.Derivative,
	}
}

func (m *rkFehlbergMethod) Forward(ivp *mned.IVP) mned.Stepper {
	return m.withStep(m.step, ivp)
}

func (m *rkFehlbergMethod) Backward(ivp *mned.IVP) mned.Stepper {
	return m.withStep(-m.step, ivp)
}

// Create an instance of the Runge-Kutta-Fehlberg method, an adaptive,
// fourth-order method, with the given parameter.
//
// This method is given by the following parameters:
// ```
// k1 = hf(t, x(t))
// k2 = hf(t+h/4, x(t)+k1/4)
// k3 = hf(t+3h/8, x(t)+3k1/32+9k2/32)
// k4 = hf(t+12h/13, x(t)+1932k1/2197-7200k2/2197+7296k3/2197)
// k5 = hf(t+h, x(t)+439k1/216-8k2+3680k3/513-845k4/4104)
// k6 = hf(t+h/2, x(t)-8k1/27+2k2-3544k3/2656+1859k4/4104-11k5/40)
// e = k1/360-128k3/4275-2197k4/75240+k5/50+2k6/55
// x(t+h) ~ x(t)+25k1/216+1408k3/2565+2197k4/4104-k5/5
// ```
// If `|e|` is greater than `tol`, we let `q:=(h*tol/(2*|e|))^(1/4)`, saturated
// into `[0.1,4]`, and redo the operations with `h*q`, saturating over `hmax`
// and stopping below `hmin`.
func RKFehlberg(
	tolerance float64, step float64, hmin float64, hmax float64,
) mned.Method {
	return &rkFehlbergMethod{
		step: step,
		tol:  tolerance,
		hmin: hmin,
		hmax: hmax,
	}
}
