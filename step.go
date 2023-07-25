// Stepping routines and method interface.
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

package mned

import (
	"math"
	"sort"
)

// A Stepper is an iterator over points of the solution of an IVP. ODE solver
// methods generally produce steppers that go from the initial values of the
// problems to lower or greater values.
type Stepper interface {
	// Go to the next point and return a pointer to it. If the second value
	// is false, there are no more points in the direction the Stepper is
	// following and the value of Point should not be trusted.
	//
	// The point returned might change after the next call to next() on the
	// same Stepper.
	Next() (*Point, bool)
}

// Starting in the point init, step until the domain of the function ends or
// one of the events says that the stepping should stop, calculating the
// intermediante points with the given interpolator and, if `callback` is not
// nil, calling the given callback with each new point. If ok is false, the
// iteration ended because the stepper returned false; otherwise index contains
// the index of the event that made the iteration stop.
func StepUntil(
	s Stepper,
	init Point,
	interp Interpolator,
	callback func(*Point),
	evs ...Event,
) (index int, ok bool) {
	var prev Point = init
	vs := make([]float64, len(evs))
	for i := 0; i < len(evs); i++ {
		vs[i] = evs[i].Cross(&prev)
	}
	for {
		point, ok := s.Next()
		if !ok {
			return 0, false
		}
		actions := make(requiredActions, 0)
		for i := 0; i < len(evs); i++ {
			newVal := evs[i].Cross(point)
			if differentSign(vs[i], newVal) {
				actions = append(actions, requiredAction{
					index: i,
					point: evs[i].FindPoint(
						interp, &prev, point,
					)})
			}
			vs[i] = newVal
		}

		sort.Sort(actions)
		size := len(actions)
		if point.Time < prev.Time {
			for i := 0; i < size/2; i++ {
				actions[i], actions[size-i] =
					actions[size-i], actions[i]
			}
		}
		for _, a := range actions {
			if !evs[a.index].Action(&a.point) {
				return a.index, true
			}
		}

		callback(point)
		prev = point.Clone()
	}
}

// A ConfigurableStepper is a Stepper which allows specifying the delta of the
// next step; that is, the difference between the time of a point and that of
// the previous point.
type ConfigurableStepper interface {
	Stepper
	// Like Stepper.Next but the time of the next point is specified to be
	// the time of the latest point plus the value given.
	NextStep(float64) (*Point, bool)
}

// A Method is a way of approximating points of the solution to an IVP. It works
// as a factory of Steppers that step forward and backward from the initial
// values of the problem.
type Method interface {
	// Make a stepper over the solution of the given IVP that starts at the
	// initial values of the problem and goes to ever increasing values for
	// the independent variable.
	Forward(*IVP) Stepper
	// Make a stepper over the solution of the given IVP that starts at the
	// initial values of the problem and goes to ever decreasing values for
	// the independent variable.
	Backward(*IVP) Stepper
}

type fixedStepper struct {
	state  IVP
	step   float64
	method func(*IVP, float64) bool
}

func (s *fixedStepper) NextStep(h float64) (*Point, bool) {
	ok := s.method(&s.state, h)
	if !ok {
		return nil, false
	}
	s.state.Start.Time += h
	return &s.state.Start, true
}

func (s *fixedStepper) Next() (*Point, bool) {
	return s.NextStep(s.step)
}

// Create a fixed step method, one where all steps increment or decrement the
// time by the same amount, with the added condition that it shouldn't depend
// on extra state. The Step is the fixed increment/decrement of time, which
// must be positive. The steppers returned by this method are
// `ConfigurableStepper`s.
//
// The Method is a function that takes an IVP and a step and stores in
// IVP.Start.Value the value for IVP.Start.Time+step, while leaving
// IVP.Start.Time untouched. The return value of `method` is true if the
// operation completed successfully or false if there was a problem, such as
// going out of the domain of the IVP.Derivative, in which case the resulting
// state of the IVP is unspecified.
type FixedStepMethod struct {
	Step   float64
	Method func(*IVP, float64) bool
}

func (m FixedStepMethod) Forward(ivp *IVP) Stepper {
	return &fixedStepper{
		state:  ivp.Clone(),
		step:   m.Step,
		method: m.Method}
}

func (m FixedStepMethod) Backward(ivp *IVP) Stepper {
	return &fixedStepper{
		state:  ivp.Clone(),
		step:   -m.Step,
		method: m.Method}
}

type adaptiveStepper struct {
	state     IVP
	step      float64
	hmin      float64
	hmax      float64
	tolerance float64
	method    func(*IVP, float64) bool
	order     uint8
}

// Return the new step that should be taken by an AdaptiveStepMethod of the
// given order and with a given tolerance took a step of a given size yielding
// the given amount of error, with max as the maximum size of a step.
func AdaptiveUpdateStep(
	step float64, tol float64, err float64, max float64, order uint8,
) float64 {
	var q float64
	adjust := tol * math.Abs(step) / (2 * err)
	if order == 1 { // Fast path
		q = adjust
	} else {
		q = math.Pow(adjust, 1/float64(order))
	}
	switch {
	case q < 0.1:
		q = 0.1
	case q > 4:
		q = 4
	}
	newStep := q * step
	if math.Abs(newStep) > max {
		if newStep > 0 {
			newStep = max
		} else {
			newStep = -max
		}
	}
	return newStep
}

func (s *adaptiveStepper) adaptStep(err float64) {
	s.step = AdaptiveUpdateStep(s.step, s.tolerance, err, s.hmax, s.order)
}

func (s *adaptiveStepper) getSteps(full []float64, half []float64) bool {
	orig := s.state.Start.Value
	copy(full, orig)
	s.state.Start.Value = full
	if !s.method(&s.state, s.step) {
		return false
	}
	copy(half, orig)
	s.state.Start.Value = half
	if !s.method(&s.state, s.step/2) {
		return false
	}
	if !s.method(&s.state, s.step/2) {
		return false
	}
	s.state.Start.Value = orig
	return true
}

func (s *adaptiveStepper) Next() (*Point, bool) {
	size := len(s.state.Start.Value)
	full := make([]float64, size)
	half := make([]float64, size)

	for s.step >= s.hmin {
		for !s.getSteps(full, half) {
			if s.step *= 0.5; s.step < s.hmin {
				return nil, false
			}
		}

		orderm1 := 1 / (math.Ldexp(1, int(s.order)) - 1)
		multiplier := 1 + orderm1
		sqerr := float64(0)
		for i, v := range full {
			diff := v - half[i]
			sqerr += diff * diff
		}
		err := math.Sqrt(sqerr) * multiplier
		tol := s.tolerance * math.Abs(s.step)
		if err > tol {
			s.adaptStep(err)
			continue
		}

		for i, v := range full {
			s.state.Start.Value[i] =
				multiplier*half[i] - orderm1*v
		}

		s.state.Start.Time += s.step
		s.adaptStep(err)
		return &s.state.Start, true
	}

	return nil, false // s.step < s.hmin
}

// An AdaptiveStepMethod is a kind of adaptive method (that is, one that changes
// the step over time for greater precision and performance) that can be adapted
// from a fixed-step method.
//
// Step, Hmin, and Hmax are the initial, minimum, and maximum step. If the
// step must go lower than Hmin to satisfy the tolerance, the stepper ends;
// if it could get greater than Hmax, it saturates to Hmax. The Tolerance is
// such that the Euclidean norm of the estimated error of any step of size `h`
// must be lower than `h*Tolerance`.
//
// Method works the same way as in `FixedStepMethod` except that its error for a
// given interval of time must be on the order of `O(h^Order)`, with the
// variable `h` being the step. Order must be positive.
//
// The algorithm used is the following, where `(t,x(t))` are the initial values,
// `y(h,t,x(t))` is the result of applying the inner method with step `h` from
// `(t,x(t))`, `k` is the order of the method and `h` is the step.
// 1. First, `F:=y(h,t,x(t))` and `H:=y(h/2,t+h/2,y(h/2,t,x(t)))` are obtained,
//    and the error is estimated as `e:=|(1+1/(2^k-1))*(H-F)|`.
// 2. If that error is lower than `|h|*Tolerance`, the step is accepted.
//    Otherwise we set `q:=(|h|*Tolerance)/(2*e)`, saturate so that `q` is in
//    `[0.1, 4]`, set `h:=q*h`, check bounds on `h` and return to the start.
// 3. The next value is taken to be `(t+h, (2^k*H - H)/(2^k - 1))`, and the
//    next step size is calculated by changing `h` as in step 2.
//
// The steppers are *not* `ConfigurableStepper`s.
type AdaptiveStepMethod struct {
	Step      float64
	Hmin      float64
	Hmax      float64
	Tolerance float64
	Method    func(*IVP, float64) bool
	Order     uint8
}

func (m *AdaptiveStepMethod) toStepper(ivp *IVP, step float64) Stepper {
	return &adaptiveStepper{
		state:     ivp.Clone(),
		step:      step,
		hmin:      m.Hmin,
		hmax:      m.Hmax,
		tolerance: m.Tolerance,
		method:    m.Method,
		order:     m.Order,
	}
}

func (m *AdaptiveStepMethod) Forward(ivp *IVP) Stepper {
	return m.toStepper(ivp, m.Step)
}

func (m *AdaptiveStepMethod) Backward(ivp *IVP) Stepper {
	return m.toStepper(ivp, -m.Step)
}
