// Correction routines for optimal step estimation.
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
	"math"
	"github.com/JwanMan/mned"
	"github.com/JwanMan/mned/method"
)

type pcStepper struct {
	steps   [5]Step
	problem mned.IVP
	h       float64
	init    uint8
}

func (s *pcStepper) NextStep(h float64) (*mned.Point, bool) {
	var ok bool
	if s.h == 0 {
		return nil, false
	}
	if s.init != 0 {
		deriv, ok := method.RK4Step(&s.problem, h)
		if !ok {
			return nil, false
		}
		s.problem.Start.Time += h
		s.steps[s.init+1].Deriv = deriv
		s.steps[s.init].Value = make([]float64, len(deriv))
		copy(s.steps[s.init].Value, s.problem.Start.Value)
		return &s.problem.Start, true
	}

	s.steps[1].Deriv, ok = s.problem.Derivative(s.problem.Start)
	if !ok {
		s.h = 0
		return nil, false
	}
	s.steps[0].Value = adamsBashfordStep(s.steps[1:], h)
	s.problem.Start.Time += h
	s.problem.Start.Value = s.steps[0].Value
	s.steps[0].Deriv, ok = s.problem.Derivative(s.problem.Start)
	if !ok {
		s.problem.Start.Time -= h
		s.problem.Start.Value = s.steps[1].Value
		return nil, false
	}
	s.steps[0].Value = adamsMoultonStep(s.steps[:], h)
	s.problem.Start.Value = s.steps[0].Value
	Shift(s.steps[:])
	s.steps[0].Value = nil
	s.steps[0].Deriv = nil
	s.steps[1].Deriv = nil
	return &s.problem.Start, true
}

func (s *pcStepper) Next() (*mned.Point, bool) {
	return s.NextStep(s.h)
}

type pcMethod float64

func (m pcMethod) withStep(h float64, ivp *mned.IVP) mned.Stepper {
	result := pcStepper{
		problem: ivp.Clone(),
		h:       h,
		init:    3,
	}
	result.steps[4].Value = make([]float64, len(ivp.Start.Value))
	copy(result.steps[4].Value, ivp.Start.Value)
	return &result
}

func (m pcMethod) Forward(ivp *mned.IVP) mned.Stepper {
	return m.withStep(float64(m), ivp)
}

func (m pcMethod) Backward(ivp *mned.IVP) mned.Stepper {
	return m.withStep(-float64(m), ivp)
}

func PredictorCorrector(step float64) mned.Method {
	return pcMethod(step)
}

type adaptivePCStepper struct {
	steps     [5]Step
	lasttime  float64
	deriv     func(mned.Point) ([]float64, bool)
	h         float64
	hmin      float64
	hmax      float64
	tolerance float64
	pos       uint8
}

func (s *adaptivePCStepper) tryStep() (errNorm float64, ok bool) {
	s.steps[1].Deriv, ok = s.deriv(mned.Point{
		Time:  s.lasttime,
		Value: s.steps[1].Value,
	})
	if !ok {
		return 0, false
	}
	s.steps[0].Value = adamsBashfordStep(s.steps[1:], s.h)

	s.steps[0].Deriv, ok = s.deriv(mned.Point{
		Time:  s.lasttime + s.h,
		Value: s.steps[0].Value,
	})
	if !ok {
		return 0, false
	}
	bash := s.steps[0].Value
	moulton := adamsMoultonStep(s.steps[:], s.h)
	s.steps[0].Value = moulton
	s.steps[0].Deriv = nil

	sqerr := float64(0)
	for i, v := range s.steps[0].Value {
		diff := v - bash[i]
		sqerr += diff * diff
	}
	err := (math.Sqrt(sqerr) * 19) / 270
	return err, true
}

func (s *adaptivePCStepper) adjustStep(err float64) {
	q := math.Pow((s.tolerance*math.Abs(s.h))/(2*err), 0.25)
	if q < 0.1 {
		q = 0.1
	}
	if q > 4 {
		q = 4
	}
	s.h *= q
	if math.Abs(s.h) > s.hmax {
		if s.h > 0 {
			s.h = s.hmax
		} else {
			s.h = -s.hmax
		}
	}
}

func (s *adaptivePCStepper) init() bool {
	var ok bool
	size := len(s.steps[4].Value)
	ivp := mned.IVP{Derivative: s.deriv}
	for math.Abs(s.h) >= s.hmin {
		ivp.Start.Time = s.lasttime
		ivp.Start.Value = make([]float64, size)
		copy(ivp.Start.Value, s.steps[4].Value)

		for i := 4; i > 1; i-- {
			s.steps[i].Deriv, ok = method.RK4Step(&ivp, s.h)
			if !ok {
				return false
			}
			ivp.Start.Time += s.h
			s.steps[i-1].Value = make([]float64, size)
			copy(s.steps[i-1].Value, ivp.Start.Value)
		}
		s.lasttime = ivp.Start.Time
		errn, ok := s.tryStep()
		s.lasttime += s.h
		if !ok {
			s.lasttime -= 4 * s.h
			return false
		}
		tol := s.tolerance * s.h
		if errn > tol {
			s.lasttime -= 4 * s.h
			s.adjustStep(errn)
			continue
		}
		s.pos = 4
		return true
	}
	return false
}

func (s *adaptivePCStepper) tryInit() bool {
	for math.Abs(s.h) >= s.hmin {
		if s.init() {
			Shift(s.steps[:])
			return true
		}
		for i := 0; i < 4; i++ {
			s.steps[i].Value = nil
			s.steps[i].Deriv = nil
		}
		s.steps[4].Deriv = nil
		s.h /= 2
	}
	return false
}

func (s *adaptivePCStepper) Next() (*mned.Point, bool) {
	var result mned.Point

	if s.pos == 5 {
		if !s.tryInit() {
			return nil, false
		}
	}
	for s.h >= s.hmin {
		if s.pos > 0 {
			result.Time = s.lasttime - s.h*float64(s.pos-1)
			result.Value = s.steps[s.pos].Value
			s.pos--
			return &result, true
		}

		errn, ok := s.tryStep()
		if !ok {
			s.h /= 2
			continue
		}
		tol := s.tolerance * s.h
		if errn > tol {
			s.steps[4] = s.steps[1]
			s.adjustStep(errn)
			if !s.tryInit() {
				return nil, false
			}
			continue
		}
		s.lasttime += s.h
		Shift(s.steps[:])
		result.Time = s.lasttime
		result.Value = s.steps[1].Value
		if 10*errn < tol {
			s.steps[4] = s.steps[1]
			s.adjustStep(errn)
			s.pos = 5
		}
		return &result, true
	}
	return nil, false
}

type adaptivePCMethod struct {
	h         float64
	hmax      float64
	hmin      float64
	tolerance float64
}

func (m *adaptivePCMethod) withStep(h float64, ivp *mned.IVP) mned.Stepper {
	result := adaptivePCStepper{
		h:         h,
		hmax:      m.hmax,
		hmin:      m.hmin,
		tolerance: m.tolerance,
		lasttime:  ivp.Start.Time,
		deriv:     ivp.Derivative,
		pos:       5,
	}
	result.steps[4].Value = make([]float64, len(ivp.Start.Value))
	copy(result.steps[4].Value, ivp.Start.Value)
	return &result
}

func (m *adaptivePCMethod) Forward(ivp *mned.IVP) mned.Stepper {
	return m.withStep(m.h, ivp)
}

func (m *adaptivePCMethod) Backward(ivp *mned.IVP) mned.Stepper {
	return m.withStep(-m.h, ivp)
}

func AdaptivePredictorCorrector(
	tolerance, step, hmin, hmax float64,
) mned.Method {
	return &adaptivePCMethod{
		h:         step,
		hmax:      hmax,
		hmin:      hmin,
		tolerance: tolerance,
	}
}
