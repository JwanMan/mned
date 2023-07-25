// Implicit and trapezium stepping methods.
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

import "github.com/JwanMan/mned"

type YDerivative func(float64, float64) float64

type implicitEulerStepper struct {
	f         func(mned.Point) ([]float64, bool)
	dfy       YDerivative
	step      float64
	tolerance float64
	current   mned.Point
}

func (s *implicitEulerStepper) NextStep(h float64) (*mned.Point, bool) {
	inner := s.f
	time := s.current.Time + s.step
	value := s.current.Value[0]
	deriv := s.dfy
	fn := func(w float64) (float64, bool) {
		v, oki := inner(mned.Point{Time: time, Value: []float64{w}})
		if oki {
			return w - value - h*v[0], true
		} else {
			return 0, false
		}
	}
	result, ok := NewtonRoot1D(
		s.current.Value[0], s.tolerance, fn,
		func(w float64) (float64, bool) {
			return 1 - h*deriv(time, w), true
		})
	if !ok {
		return nil, false
	}
	s.current.Value[0] = result
	s.current.Time += s.step
	return &s.current, true
}

func (s *implicitEulerStepper) Next() (*mned.Point, bool) {
	return s.NextStep(s.step)
}

type implicitEulerMethod struct {
	dfy       YDerivative
	h         float64
	tolerance float64
}

func (m *implicitEulerMethod) withStep(h float64, ivp *mned.IVP) mned.Stepper {
	return &implicitEulerStepper{
		f:         ivp.Derivative,
		current:   ivp.Start.Clone(),
		step:      h,
		dfy:       m.dfy,
		tolerance: m.tolerance,
	}
}

func (m *implicitEulerMethod) Forward(ivp *mned.IVP) mned.Stepper {
	return m.withStep(m.h, ivp)
}

func (m *implicitEulerMethod) Backward(ivp *mned.IVP) mned.Stepper {
	return m.withStep(-m.h, ivp)
}

func ImplicitEuler(step, tolerance float64, df YDerivative) mned.Method {
	return &implicitEulerMethod{h: step, tolerance: tolerance, dfy: df}
}

type trapeziumStepper struct {
	f         func(mned.Point) ([]float64, bool)
	dfy       YDerivative
	step      float64
	tolerance float64
	current   mned.Point
	prevf     float64
}

func (s *trapeziumStepper) NextStep(h float64) (*mned.Point, bool) {
	inner := s.f
	time := s.current.Time + s.step
	pv, ok := s.f(mned.Point{
		Time:  s.current.Time,
		Value: []float64{s.current.Value[0]},
	})
	if !ok {
		return nil, false
	}
	value := s.current.Value[0] + h/2*pv[0]
	fn := func(w float64) (float64, bool) {
		v, oki := inner(mned.Point{Time: time, Value: []float64{w}})
		if oki {
			return w - value - h/2*v[0], true
		} else {
			return 0, false
		}
	}
	deriv := s.dfy
	result, ok := NewtonRoot1D(
		s.current.Value[0], s.tolerance, fn,
		func(w float64) (float64, bool) {
			return 1 - h/2*deriv(time, w), true
		})
	if !ok {
		return nil, false
	}
	s.current.Value[0] = result
	s.current.Time += s.step
	return &s.current, true
}

func (s *trapeziumStepper) Next() (*mned.Point, bool) {
	return s.NextStep(s.step)
}

type trapeziumMethod struct {
	dfy       YDerivative
	h         float64
	tolerance float64
}

func (m *trapeziumMethod) withStep(h float64, ivp *mned.IVP) mned.Stepper {
	return &trapeziumStepper{
		f:         ivp.Derivative,
		current:   ivp.Start.Clone(),
		step:      h,
		dfy:       m.dfy,
		tolerance: m.tolerance,
	}
}

func (m *trapeziumMethod) Forward(ivp *mned.IVP) mned.Stepper {
	return m.withStep(m.h, ivp)
}

func (m *trapeziumMethod) Backward(ivp *mned.IVP) mned.Stepper {
	return m.withStep(-m.h, ivp)
}

func Trapezium(step, tolerance float64, df YDerivative) mned.Method {
	return &trapeziumMethod{h: step, tolerance: tolerance, dfy: df}
}
