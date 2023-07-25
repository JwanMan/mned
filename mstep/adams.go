// Adams stepping methods.
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

package mstep

import (
	"github.com/JwanMan/mned"
	"github.com/JwanMan/mned/method"
)

func adamsBashfordStep(steps []Step, h float64) []float64 {
	return Implicit(steps, []float64{1}, []float64{
		float64(55) / 24, float64(-59) / 24,
		float64(37) / 24, float64(-9) / 24,
	}, h)
}

func adamsMoultonStep(steps []Step, h float64) []float64 {
	return Implicit(steps, []float64{0, 1}, []float64{
		float64(9) / 24, float64(19) / 24,
		float64(-5) / 24, float64(1) / 24,
	}, h)
}

type adamsBashfordStepper struct {
	steps   [4]Step
	h       float64
	problem mned.IVP
	init    uint8
}

func (s *adamsBashfordStepper) NextStep(h float64) (*mned.Point, bool) {
	var ok bool
	if s.h == 0 { // Error indicator
		return nil, false
	}
	if s.init != 0 {
		deriv, ok := method.RK4Step(&s.problem, h)
		if !ok {
			return nil, false
		}
		s.problem.Start.Time += h
		s.steps[s.init].Deriv = deriv
		s.init--
		s.steps[s.init].Value = make([]float64, len(deriv))
		copy(s.steps[s.init].Value, s.problem.Start.Value)
		return &s.problem.Start, true
	}

	s.steps[0].Deriv, ok = s.problem.Derivative(s.problem.Start)
	if !ok {
		s.h = 0
		return nil, false
	}
	s.problem.Start.Value = adamsBashfordStep(s.steps[:], h)
	s.problem.Start.Time += h
	Shift(s.steps[:])
	s.steps[0].Value = s.problem.Start.Value
	return &s.problem.Start, true
}

func (s *adamsBashfordStepper) Next() (*mned.Point, bool) {
	return s.NextStep(s.h)
}

type adamsBashfordMethod float64

func (m adamsBashfordMethod) withStep(h float64, ivp *mned.IVP) mned.Stepper {
	result := adamsBashfordStepper{
		h:       h,
		problem: ivp.Clone(),
		init:    3,
	}
	result.steps[3].Value = make([]float64, len(ivp.Start.Value))
	copy(result.steps[3].Value, ivp.Start.Value)
	return &result
}

func (m adamsBashfordMethod) Forward(ivp *mned.IVP) mned.Stepper {
	return m.withStep(float64(m), ivp)
}

func (m adamsBashfordMethod) Backward(ivp *mned.IVP) mned.Stepper {
	return m.withStep(-float64(m), ivp)
}

func AdamsBashford(step float64) mned.Method {
	return adamsBashfordMethod(step)
}
