// Core data type definitions.
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

/*
This package implements several methods to compute approximations of
solutions of initial value problems (IVPs) with first-order ODEs. These problems
have the form
```
x'(t) = f(t, x(t)),
x(t0) = x0,
```
where `t0 : R` and `x0 : R^n` are the **initial values**,
`f : Ω ⊆ R x R^n -> R^n` is the **derivative function**, and `x : I ⊆ R -> R^n`
is the unknown of the equation. Here `I` is an open interval of `R` containing
`t0` and `Ω` is an open subset of `R x R^n` containing `(t0, x0)`.

We refer to `t` as the **independent variable** and to `x(t)` as the `dependent
variable`. We know that `x` is uniquely defined in a neighbourhood of `t0`.
Given a point `t` in such a neighbourhood, we say that `(t,x(t))` is a point of
the solution of the initial value problem.

The results calculated will be approximations, and judgement is required to get
approximations suitable for a particular problem; nevertheless, for the purpose
of explanation, we will use the same terminology for the approximations than for
the actual functions.
*/

// A Point is an element of the solution, given by the values of the independent
// and dependent variables.
type Point struct {
	Time  float64   // The independent variable.
	Value []float64 // The dependent variable.
}

// Deep-copy a point.
func (p *Point) Clone() Point {
	value := make([]float64, len(p.Value))
	copy(value, p.Value)
	return Point{Time: p.Time, Value: value}
}

// An IVP is an initial value problem given by a first-order ODE. It's given by
// the Derivative function, which is a pure function, and the initial values.
//
// The second return value of Derivative indicates if the Point was in the
// domain of the function. If it wasn't the first return value should not be
// used. Ideally, the domain should be restricted to the connected component of
// the initial value, as otherwise a solving method could "jump" to another
// component and the results from there on would be invalid.
type IVP struct {
	Derivative func(Point) ([]float64, bool)
	Start      Point
}

func (ivp *IVP) Clone() IVP {
	return IVP{Derivative: ivp.Derivative, Start: ivp.Start.Clone()}
}
