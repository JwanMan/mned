// Caching of solution points.
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

import "sort"

// A CacheSolver is like a dynamic version of a DenseSolution. It stores
// points of the solution of a problem in a conceptually maximal interval that
// is evaluated lazily; that is, it starts with the solutions in a range that
// only includes the initial point and, when a point is asked for a Time value
// outside the range, the range is expanded to cover that point.
type CacheSolver struct {
	backward []Point // Invariant: Decreasing times, times before forward's.
	forward  []Point // Invariant: Non-empty, increasing times.
	backStep Stepper // Invariant: Decreasing times from end of backward.
	forStep  Stepper // Invariant: Increasing times from end of forward.
	evs      []Event
	interp   Interpolator
}

// Create a CacheSolver to solve the given IVP with the given solving Method,
// with the given interpolator to calculate points between steps and taking into
// account the events provided.
func CacheSolve(
	m Method, ivp *IVP, interp Interpolator, evs ...Event,
) CacheSolver {
	return CacheSolver{
		forward:  []Point{ivp.Start.Clone()},
		backward: []Point{},
		forStep:  m.Forward(ivp),
		backStep: m.Backward(ivp),
		evs:      evs,
		interp:   interp,
	}
}

// Get the time of the stored solution point with the lowest time.
func (c *CacheSolver) Start() float64 {
	if len(c.backward) > 0 {
		return c.backward[len(c.backward)-1].Time
	} else {
		return c.forward[0].Time
	}
}

// Get the time of the stored solution point with the greatest time.
func (c *CacheSolver) End() float64 {
	return c.forward[len(c.forward)-1].Time
}

func (c *CacheSolver) getActions(
	current *Point, next *Point, vals []float64,
) requiredActions {
	var actions requiredActions = make([]requiredAction, 0)
	for i := 0; i < len(vals); i++ {
		newVal := c.evs[i].Cross(next)
		if differentSign(vals[i], newVal) {
			point := c.evs[i].FindPoint(c.interp, current, next)
			actions = append(
				actions,
				requiredAction{index: i, point: point},
			)
		}
		vals[i] = newVal
	}
	sort.Sort(actions)
	return actions
}

func (c *CacheSolver) expandBackward(t float64) bool {
	if c.backStep == nil {
		return false
	}

	var current Point
	if len(c.backward) == 0 {
		current = c.forward[0]
	} else {
		current = c.backward[len(c.backward)-1]
	}

	vals := make([]float64, len(c.evs))
	for i := 0; i < len(vals); i++ {
		vals[i] = c.evs[i].Cross(&current)
	}

	for t < current.Time {
		next, ok := c.backStep.Next()
		if !ok {
			// TODO Consider bisection for reasonable last point
			c.backStep = nil
			return false
		}
		actions := c.getActions(&current, next, vals)
		for i := len(actions) - 1; i >= 0; i-- {
			if !c.evs[actions[i].index].Action(&actions[i].point) {
				return false
			}
		}
		current = next.Clone()
		c.backward = append(c.backward, current)
	}
	return true
}

func (c *CacheSolver) expandForward(t float64) bool {
	if c.forStep == nil {
		return false
	}

	current := c.forward[len(c.forward)-1]
	vals := make([]float64, len(c.evs))
	for i := 0; i < len(vals); i++ {
		vals[i] = c.evs[i].Cross(&current)
	}

	for t < current.Time {
		next, ok := c.forStep.Next()
		if !ok {
			// TODO Consider bisection for reasonable last point
			c.forStep = nil
			return false
		}
		actions := c.getActions(&current, next, vals)
		for i := 0; i < len(actions); i++ {
			if !c.evs[actions[i].index].Action(&actions[i].point) {
				return false
			}
		}
		current = next.Clone()
		c.forward = append(c.forward, current)
	}
	return true
}

// Get the value of the solution for a given time.
//
// If the value is outside of the bounds of what's stored, values are added
// until reaching the given time. If ok is false, the value is out of bounds for
// the problem or some event and points have only been added as far as it was
// possible.
func (c *CacheSolver) Get(t float64) (x []float64, ok bool) {
	if t < c.Start() {
		if ok := c.expandBackward(t); !ok {
			return nil, false
		}
	}
	if t > c.End() {
		if ok := c.expandForward(t); !ok {
			return nil, false
		}
	}
	if t >= c.forward[0].Time {
		i := 1
		for t > c.forward[i].Time {
			i++
		}
		return c.interp.FindValue(
			&c.forward[i-1], &c.forward[i], t,
		), true
	}

	i := 0
	for t > c.backward[i].Time {
		i++
	}
	if i == 0 {
		return c.interp.FindValue(
			&c.backward[0], &c.forward[0], t,
		), true
	}
	return c.interp.FindValue(&c.backward[i], &c.backward[i-1], t), true
}

// Get the value at `c.Start() - h` as in `c.Get` except that, if the underlying
// backward Stepper is a ConfigurableStepper, its method `NextStep(-h)` is used.
// The result value should *NOT* be mutated.
//
// This is most useful with fixed step methods to specify the step.
func (c *CacheSolver) StepBackward(h float64) ([]float64, bool) {
	if cs, ok := c.backStep.(ConfigurableStepper); ok {
		if point, ok := cs.NextStep(-h); ok {
			if h > 0 {
				c.backward = append(c.backward, point.Clone())
			}
			return point.Value, true
		} else {
			return nil, false
		}
	}
	return c.Get(c.Start() - h)
}

// Get the value at `c.End() + h` as in `c.Get` except that, if the underlying
// forward Stepper is a ConfigurableStepper, its method `NextStep(h)` is used.
// The result value should *NOT* be mutated.
//
// This is most useful with fixed step methods to specify the step.
func (c *CacheSolver) StepForward(h float64) ([]float64, bool) {
	if cs, ok := c.forStep.(ConfigurableStepper); ok {
		if point, ok := cs.NextStep(h); ok {
			if h > 0 {
				c.forward = append(c.forward, point.Clone())
			}
			return point.Value, true
		} else {
			return nil, false
		}
	}
	return c.Get(c.End() + h)
}

// Get the points calculated from the initial values to earlier times. The
// points **MUST NOT** be modified.
func (c *CacheSolver) BackwardPoints() []Point {
	return c.backward
}

// Get the points calculated from the initial values to later times. The points
// **MUST NOT** be modified.
func (c *CacheSolver) ForwardPoints() []Point {
	return c.forward
}

// Rearrange the stored solution points in a result such that the i-th point has
// time `result[0][i]` and value `(result[1][i],...,result[size][i])`, where
// `size` is the number of dimensions of a value. This is useful for plotting.
func (c *CacheSolver) PointCoords() [][]float64 {
	size := len(c.forward[0].Value)
	result := make([][]float64, size+1)
	backs := len(c.backward)
	elems := backs + len(c.forward)

	for i := 0; i <= size; i++ {
		result[i] = make([]float64, elems)
	}
	for i := 0; i < backs; i++ {
		result[0][i] = c.backward[backs-i-1].Time
		for j := 0; j < size; j++ {
			result[j+1][i] = c.backward[backs-i-1].Value[j]
		}
	}
	for i := backs; i < elems; i++ {
		result[0][i] = c.forward[i-backs].Time
		for j := 0; j < size; j++ {
			result[j+1][i] = c.forward[i-backs].Value[j]
		}
	}
	return result
}

// Force the generation of values to the left until the domain ends or an event
// tells the stepping to stop.
func (c *CacheSolver) StepToBeginning() {
	var init Point
	if len(c.backward) > 0 {
		init = c.backward[len(c.backward)-1]
	} else {
		init = c.forward[0]
	}
	StepUntil(
		c.backStep,
		init,
		c.interp,
		func(point *Point) {
			c.backward = append(c.backward, point.Clone())
		},
		c.evs...,
	)
}

// Force the generation of values to the right until the domain ends or an event
// tells the stepping to stop.
func (c *CacheSolver) StepToEnd() {
	StepUntil(
		c.forStep,
		c.forward[len(c.forward)-1],
		c.interp,
		func(point *Point) {
			c.forward = append(c.forward, point.Clone())
		},
		c.evs...,
	)
}
