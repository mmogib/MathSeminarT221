### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ a563c4b9-01d0-4411-97ae-8019414330b1
begin
	using Random
	using Statistics
	using LinearAlgebra
	using JuMP, Ipopt 
	using HiGHS
	using PlutoUI
	using Plots
	using StatsPlots
	using LaTeXStrings
	using HTTP, JSON3
end

# ╔═╡ c154e809-bee8-444f-921f-32100509519e
TableOfContents(title="MATH Seminar, Monday Oct. 31")

# ╔═╡ 07d5b126-4150-11ed-24e1-eb98b921c8c6
md"""
# Optimization Algebraic Modeling
## Abstract
Many researchers at a point of time have to solve an optimization problem. These problems range from simple unconstrained function to a complex model of many variables and constraints. In this talk, we introduce the concept of algebraic modeling in optimization which is used to translate between the algebraic form of the modeler and the standard form of the solver. It aims to allow researchers to express their optimization problems in a natural, algebraic form similar to the original mathematical expressions. Examples from different branches of optimization will be given for illustration.
"""

# ╔═╡ d716ec58-253d-4267-9497-7497fe639c6a
md"""
__Questions__ Before this seminar,
1. Have you heard of __*algebraic modeling in optimization*__ and __*algebraic modeling languages*__? 
2. Have you ever used __*optimization techniques*__ to solve problems in your work?
3. Have you used __*Julia*__ programming language?
4. Have you used __*JuMP*__ modeling language?

__Objective of this talk__

- Introduce __*algebraic modeling in optimization*__ and __*algebraic modeling languages*__.
- Give general guidelines on how any researcher can use state-of-the-art solvers with minimal effort.
- Introduce __*JuMP*__ modeling language with __examples__.
- Discuss some __challenges__ that are suitable for an MS/PhD thesis in Mathematict/Computer Science.

"""

# ╔═╡ bcfd72c1-c3f5-4cc6-96dd-b9d292d78c7f
md"""
## Introduction
A general __optimization problem__ is one of the form
```math
\begin{array}{ll}
\min & f(x) \\
\text{subject to}\\
 & x\in \Omega.
\end{array} \tag{OP}
```
A model in mathematical optimization consists of four key objects:
- __data__ [also called the constants of a model], 
- __variables__ (continuous, semi-continuous, binary, integer) [also called decision variables or parameters],
- __constraints__ (equalities, inequalities) [sometimes also called restrictions], and
- the __objective function__ [also called the cost function]

Depending on ``f``, ``x``, and ``\Omega``, we have different types of problems, such as 
- linear programming (__LP__) problems
- mixed integer linear programming (__MILP__) problems, 
- nonlinear programming (__NLP__) problems, and,
- mixed integer nonlinear programming (__MINLP__) problems.
- vector optimization (__VO__) problems.

"""

# ╔═╡ f111a4cf-7f17-429b-9bb0-290fe484e2d0
md"""
## Example
Here's the problem:
```math
\begin{aligned}
& \min & 12x + 20y \\
& \;\;\text{s.t.} & 6x + 8y \geq 100 \\
& & 7x + 12y \geq 120 \\
& & x \geq 0 \\
& & y \in [0, 3] \\
\end{aligned}
```

"""

# ╔═╡ b193fb26-5736-4b37-9c8c-6a44072b4a6d
md"""
To solve this problem, we need to employ a solver. One such solve is [HiGHS](https://highs.dev/). This is an Open source serial and parallel solvers for large-scale sparse linear programming (LP), mixed-integer programming (MIP), and quadratic programming (QP) models developed in C++.

To use this solver, we must __load the model__. Loading a model in HiGHS is via
- from data file (.lp, .mps)
- via data provided by another application


```
NAME          example.mps
ROWS
 N  obj     
 G  c1      
 G  c2      
COLUMNS
    x        obj                  12   c1                   6
    x        c2                   7
    y        obj                  20   c1                   8
    y        c2                   12
RHS
    rhs       c1                  100  c2                  120
BOUNDS
 LO BOUND     x                  0
 UP BOUND     y                  3
 LO BOUND     y                  0
ENDATA
```
"""

# ╔═╡ 01488a29-905c-4ea4-8fee-d0c49b81f1e0
md"""
However, we can use an algebraic modeling language to write the problem and pass it through to the solver. We use JuMP 
```julia
model1 = Model(HiGHS.Optimizer)
set_silent(model1)
@variable(model1, x >= 0)
@variable(model1, 0 <= y <= 3)
@objective(model1, Min, 12x + 20y)
@constraint(model1, c1, 6x + 8y >= 100)
@constraint(model1, c2, 7x + 12y >= 120)
optimize!(model1);
termination_status(model1)
primal_status(model1)
dual_status(model1)
obj1v = objective_value(model1);
xv = value(x)
yv = value(y)
c1v = shadow_price(c1)
c2v = shadow_price(c2)
```
"""

# ╔═╡ 111e9ba3-7a25-4a27-90f3-fe358d1a66a5
begin
	c11RHS = @bind c11RHSv NumberField(90 : 0.1 : 110, default=100)
	c12RHS = @bind c12RHSv NumberField(110 : 0.1 : 130, default=120)
	md"""
	|RHS of constraint 1| RHS of constraint 2|
	|---|---|
	| $c11RHS | $c12RHS |
	"""
end

# ╔═╡ 3327bf8e-d5d3-44f5-85c1-c1bde64abcc3
L"""
\begin{aligned}
& \min & 12x + 20y \\
& \text{subject to} \\
& & 6x + 8y \geq %$c11RHSv \\
& & 7x + 12y \geq %$c12RHSv \\
& & x \geq 0 \\
& & y \in [0, 3] \\
\end{aligned}
"""

# ╔═╡ 5b57a081-d9e7-42d3-b579-8cf579d851bf
begin 
	model1 = Model(HiGHS.Optimizer)
	set_silent(model1)
	@variable(model1, x1 >= 0)
	@variable(model1, 0 <= y1 <= 3)
	@objective(model1, Min, 12x1 + 20y1)
	@constraint(model1, c11, 6x1 + 8y1 >= c11RHSv)
	@constraint(model1, c12, 7x1 + 12y1 >= c12RHSv)
	optimize!(model1);
	t1_status = termination_status(model1)
	p1_status = primal_status(model1)
	d1_status = dual_status(model1)
	obj1v = objective_value(model1);
	x1v = value(x1)
	y1v = value(y1)
	c11v = shadow_price(c11)
	c12v = shadow_price(c12)
	md"""
	__objective value = $(objective_value(model1))__

	__termination status__ $(t1_status)
	
	__primal status__ $(p1_status)
	
	__dual status__ $(d1_status)
	
	``x`` = $x1v
	
	``y`` = $y1v

	shadow price 1 : __$c11v__ , 	right hand side = $(6*x1v+8*y1v)
	
	shadow price 2 : __$c12v__ , 	right hand side = $(7*x1v+12*y1v)
	
	
	"""
end

# ╔═╡ 0e4cdd43-3fc0-467f-b004-250fbaf53ac0
md"""
## History
- The __first__ computer __LP solver__ was in the late __1950__'s when William OrchardHays, working at the __RAND Corporation__ in consultation with __George Dantzig__, developed a system with the computer code on __punched cards__.
- In the __1960__'s, supported by oil companies, __MPS/360__ (later MPSX) was devolped. (MPS: Mathematical Programming System)
- In __1970__'s General Algebraic Modeling System [(GAMS)](https://www.gams.com/latest/docs/gams.pdf) was designed and introduced By Meeraus and Bisschop while working for the World Bank.
- In __1983__, __mp-model__ (which is a declarative modeling language for LP and MILP problems). This was replaced by [(Mosel)](https://www.fico.com/fico-xpress-optimization/docs/latest/mosel/mosel_lang/dhtml/moselreflang.html) in 2001.
- In 1985-86, at Bell Laboratories, Robert Fourer, David M. Gay, and Brian W. Kernighan developed __AMPL__.
- In 1989, [LINGO](https://www.lindo.com/index.php/products/lingo-and-optimization-modeling) appeared.
- In early __1990__s, __AIMMS__ was deveolped as the first modeling system to offer graphical user interfaces.
- In __2006__, __[CVX](#cvx2006)__, which is a package for specifying and solving convex programs written in __MATLAB__.
- In __2008__, the first version of __Pyomo__ appeared. Pyomo (Python Optimization Modeling Objects) is a collection of packages for the general purpose programming language Python and provides tools and constructs needed to build models, declare variables, sets, solver interfaces, etc.
- In 2015, __JuMP__ first prototype appeared. It is one of, if not, the youngest environment to create optimization models. It is written in Julia programming language (first version was released in 2012).

"""

# ╔═╡ 85c4eb73-c2b2-40e3-85f8-09d5639f5d37
md"""
## Modeling Systems

1. __Algebraic Modeling Languages__
2. __Non-algebraic Modeling Languages__: build models in an __object oriented way__.
3. __Integrated Modeling Environments__: A more __graphical approach__ to model building is taken by the so called integrated modeling environments. exploiting graphical user interfaces (GUI) for various tasks.
4. __Model-Programming Languages__: __built-in scripting capabilities__, but these new languages provide elements for declaring models and for describing solution algorithms in a much more integrated way. 
"""

# ╔═╡ b51da029-f797-45d9-9309-c18050474ca6
md"""
### Algebraic Modeling Languages: 
An __algebraic modeling langauge (AML)__ takes a __mathematical optimization problem__ written by a
user, converts it into __a standard form__, passes that standard form to a __solver__, waits for the solver to complete, then queries the solver for a solution and __returns the solution to the user__.

__Algebraic modeling languages__ [[wiki202208](#wiki202208)] are high-level computer programming languages for describing and solving high complexity problems for large scale mathematical computation.

Some algebraic modeling languages such as
- AIMMS
- AMPL
- GAMS
- Gekko
- MathProg
- Mosel
- OPL
- CVX
- Pyomo
- JuMP
have syntax similar to the mathematical notation of optimization problems.


Algebraic modeling languages can be viewed as a __new paradigm__ of programming. Usually programming languages are divided into three different classes:
- __Imperative Languages__: procedural languages and object oriented programming, for example C, C++, Pascal, FORTRAN, Java.
- __Functional Languages__:  every computation can be viewed as a function ``f : X \to Y`` translating an input from ``X`` to a unique output in ``Y``. For example: LISP, haskel.
- __Logic Programming Languages__:It was designed to prove mathematical theorems. Every mathematical proof can be regarded as a computation following specific rules. In the reverse, every computation can as well be regarded as a proof. For example Prolog.

"""

# ╔═╡ ab16ed26-4aa7-40d7-8013-a2ed51352624
md"""
To date, __AMPL__, __GAMS__, and similar commercial packages represent __the state of the art__ in AMLs and are widely used in both __academia__ and __industry__. These AMLs are quite efficient at what they were designed for; however, they have a number of __drawbacks__.
- They are relatively __standalone commercial systems__ and difficult to fit naturally within a __modern scientific workflow__
  - solving optimization problems within a __larger simulation__
  - solving optimization problems within an __interactive visualization__
  - constructing a __complex optimization model programmatically__ from modular components
- It is __tedious__ for __algorithm developers__ to interact with __solvers__ while they are running, both for control of the solution process and to reduce the overhead of regenerating a model when solving a sequence of related instances.
- __Modelers__ want to __create user-friendly AML extensions__ for __new problem classes__ that couple with specialized solution approaches. These are NOT designed to be extended in this way

"""

# ╔═╡ ab450756-5cfc-4df3-9bb2-ad0b9da8c948
md"""
## JuMP
[JuMP](#jump2017) is an open-source modeling language that allows users to express a wide range of
optimization problems (linear, mixed-integer, quadratic, conic-quadratic, semidefinite, and nonlinear) in a high-level, algebraic syntax. 

- JuMP is an AML which is embedded in the __Julia programming language__ and can be easily integrated within a scientific workflow.
- Provides a __performant open-source__ alternative to commercial systems.
-  JuMP is easily __extensible__.
- JuMP was first __released in 2013__ to support linear and mixed-integer optimization.
- JuMP enables modeling standard problem classes supported by the commercial and open-source AMLs. 
  - quadratic, 
  - conic-quadratic, 
  - semidefinite, and 
  - general derivative-based nonlinear optimization problems, 
- JuMP extended beyond what is typically available in an AML, either commercial or open-source.
  - callbacks for in-memory bidirectional communication with branchand-bound solvers
  - AD of user-defined nonlinear functions,
  - easy-to-develop addons for specialized problem classes such as robust optimization (RO) and Vector Optimization (VO).
"""

# ╔═╡ 46e87d72-d899-43dd-9c81-7eb01200e6e7
md"""
## Computing Derivatives for Nonlinear Models
The role of a modeling language for __nonlinear optimization__ is to 
- allow users to specify “closed-form” algebraic expressions for the objective function and constraints and 
- communicate first-order and typically second-order derivatives with the optimization solver. 

__JuMP__, like __AMPL__ and __GAMS__, uses techniques from __automatic (or algorithmic)__ differentiation to __evaluate derivatives__ of __user-defined expressions__.

### Expression Graphs and Reverse-Mode AD.
- the expression graph, which is a __directed acyclic graph__ (or, typically, a __tree__) that encodes the sequence of operations required to compute the expression as well as the dependency structure between operations.
- For example the expression ``e^{x^2+y^2}`` is represented by $(Resource("https://www.dropbox.com/s/50m0nnqhfrkz6o5/expression_graph.png?raw=1", :height=>250)) 
   - with nodes representing the input values ``x`` and ``y`` together with every __“basic” operation__ like addition and exponentiation that is performed in computing the value of the expression. 
   - Edges in the expression graph represent immediate dependencies between operations.
   - The expression graph encodes all needed information for JuMP to evaluate and compute derivatives of nonlinear expressions, and JuMP generates these objects by using __macros__ (meta-programming).

Given an expression graph object, 
- we compute the numerical value of the expression by iterating through the nodes of the graph in an order such that by the time we reach a given node to evaluate its corresponding operation, the numerical values of all its inputs (children) have already been computed. 
- We apply the chain rule in such a way that by iterating through the nodes in the reverse order (parents before children), in a single pass, we obtain the exact gradient vector __∇f(x)__. This reverse-pass algorithm, known suitably as __reverse-mode AD__

"""


# ╔═╡ c57bcac2-9cb7-4c6d-ac7f-3436d3cbcc6f
md"""
## User-defined functions and Derivative-free
In the cases where the __built-in library is insufficient__, there has historically been no user-friendly way to incorporate user-defined functions into AMLs.

A compelling application for user-defined
functions is __optimal control problems constrained by differential equations__

JuMP is the first __AML__ to not only 
- provide a simple interface for __user-defined functions__ with __user-provided (hand-coded) derivatives__, but also to 
- provide an option to __automatically differentiate user-defined functions__.

__Example__
```math
\begin{aligned}
& \max & x_1 + x_2 \\
& \text{subject to} \\
& & \sqrt{x_1^2+x_2^2} \leq 1 \\
\end{aligned}
```
"""

# ╔═╡ 43a9716e-d735-47d1-b209-77754bad82ab
starting_point = @bind start_point_v Slider(-2:0.1:2; default=0.5, show_value=true)

# ╔═╡ 78e26d05-cd83-4029-a300-b21a38f7d724
function my_squareroot(x)
	z=x# Initial starting point for Newton’s method
	while abs(z*z - x) > 1e-13
	z = z - (z*z-x)/(2z)
	end
	return z
end	

# ╔═╡ 738f6aa6-7780-41c4-89ef-948946ba81b6
begin
	model2 = Model(Ipopt.Optimizer)
	register(model2,:my_squareroot, 1, my_squareroot, autodiff=true)
	set_silent(model2)
	@variable(model2,x[1:2], start=start_point_v)
	@objective(model2,Max,sum(x))
	@NLconstraint(model2,my_squareroot(x[1]^2+x[2]^2)<=1)
	# @NLconstraint(model2,sqrt(x[1]^2+x[2]^2)<=1)
	optimize!(model2)
	t2_status = termination_status(model2)
	p2_status = primal_status(model2)
	d2_status = dual_status(model2)
	obj2v = objective_value(model2);
	xv = value.(x)
	md"""
	__objective value = $(objective_value(model2))__

	__termination status__ $(t2_status)
	
	__primal status__ $(p2_status)
	
	__dual status__ $(d2_status)
	
	``x_1`` = $(xv[1])
	
	``x_2`` = $(xv[2])
	
	"""
	
end

# ╔═╡ 197a354d-e466-4aeb-8590-b176f6fc96bd
solution_summary(model2, verbose=true)

# ╔═╡ 38573590-c70a-413c-a141-8733b0ee1a6e
md"""
__Remarks__

- The function __```register```__ 
  - registers the nonlinear function with the symbolic name __```my_squareroot```__ and passes a reference to the function defined above. 
  - The second argument __1__ indicates that the input to the function is __univariate__. 
  - The __```autodiff=true```__ option instructs JuMP to automatically __compute the derivatives__ of the user-defined function.

"""

# ╔═╡ 53357188-bb77-41aa-b7a5-3dfe7ff4739e
md"""
🔍 No mature implementation of this approach existed in Julia of __reverse-mode AD__
"""

# ╔═╡ 57855292-adf7-4654-9583-c76d241f055b
md"""
__JuMP__ uses __[ForwardDiff.jl](#fad2016)__, a standalone Julia package implementing forward-mode AD, to compute derivatives of these functions.
"""

# ╔═╡ 745f32fb-1f54-4cd2-b6ff-73df903165a3
md"""
## Extension for Robust Optimization
Robust optimization (RO) is a methodology for addressing uncertainty in optimization problems that has grown in popularity over the last decade.

We solve an __RO problem__ with respect to the worst-case realization of those uncertain parameters over their __uncertainty set__, i.e.,
```math
\begin{array}{lll}
\displaystyle \min_{x\in X} & f(x) \\
\text{subjet to} \\
& g(x,\zeta) \leq 0, & \quad \text{for all }\zeta \in U 
\end{array}
```
where 
- ``x`` are the decision variables,
- ``\zeta`` are the uncertain parameters drawn from the uncertainty set ``U``.
- ``f: X\to \mathbb{R}`` is a function of ``x`` only.
- ``g:X\times U \to \mathbb{R}^n`` is a vector-valued function of both ``x`` and ``\zeta``.

[JuMPeR](#ro2016) is an extension for JuMP that enables modeling RO problems directly by introducing the Uncertain modeling primitive for uncertain parameters.
The syntax is essentially unchanged from JuMP, except that constraints containing
only Uncertains and constants are treated distinctly from other constraints as they
are used to define the uncertainty set.
"""

# ╔═╡ c769c0ab-9f28-4939-812e-4329a8242d14
md"""
## (Challenge/Opportunity) Vector Optimization
In the optimization problem. 

```math
\begin{array}{ll}
\min & f(x) \\
\text{subject to}\\
 & x\in \Omega,
\end{array}\tag{VOP}
```

let ``f: \mathbb{R}^n\to \mathbb{R}^m`` be a given vector-valued funcation, ``\emptyset\neq\Omega\subseteq \mathbb{R}^n``. This is a __vector optimization problem (VOP)__

🔍 JuMP does support __VOP__. 

The challenges are 
- Passing multiple objectives to the solver
- Defining order
- Querying results from the solver.

"""

# ╔═╡ 76a3d98d-c700-476e-92a7-a4d5cbf686cf
md"""
## Sudoku
__Sudoku__ is a popular number puzzle. The goal is to place the digits 1,...,9 on a nine-by-nine grid, with some of the digits already filled in. Your solution must satisfy the following rules:

- The numbers 1 to 9 must appear in each 3x3 square
- The numbers 1 to 9 must appear in each row
- The numbers 1 to 9 must appear in each column
Here is a partially solved Sudoku problem:
"""

# ╔═╡ c0c3a7dd-5a54-4966-b980-c786a6861ede
begin
 	
	sdkoSelect = @bind sdkoLevel Select(["Easy", "Medium","Hard"], default="Easy")
	md"""
	
	"""
end

# ╔═╡ d6fa6332-d7c6-4283-906e-b25b66fa6067
begin
	btnSkdo = @bind sdkoBtn Button("Give me a new game")
	md"""
	$btnSkdo
	Choose level of difficulty 
	$sdkoSelect
	"""
end


# ╔═╡ 296d6c22-0ac1-4003-b93e-6d0261402e3c
html"""
<style>
table, td, th {  
  border: 1px solid #ddd;
  text-align: left;
}

table {
  border-collapse: collapse;
  width: 80%;
}

th, td {
  padding: 5px;
 text-align: center;
}
</style>
"""


# ╔═╡ 06c991d3-c550-44b9-99b1-563c81088066
md"""
Solving a __Sudoku is not an optimization problem__ with an objective; its actually a __feasibility problem__: we wish to find a feasible solution that satisfies these rules. You can think of it as __an optimization problem with an objective of 0__.

__Let's model this problem__

⭕ Let ``x_{i,j,k}`` be a __binary variable__ for every possible number in each possible cell. [``x_{ijk} = 1`` if and only if cell ``(i,j)`` has number ``k``, where ``i`` is the row and ``j`` is the column.]
```julia
@variable(sudoku, x[i = 1:9, j = 1:9, k = 1:9], Bin);
```

⭕ Constraints

__Constraint 1__:there can be only one number per cell.
```julia
for i in 1:9  # For each row
    for j in 1:9  # and each column
        # Sum across all the possible digits. One and only one of the digits
        # can be in this cell, so the sum must be equal to one.
        @constraint(sudoku, sum(x[i, j, k] for k in 1:9) == 1)
    end
end
```

__Constraint 2__: for the rows and the columns
```julia
for ind in 1:9  # Each row, OR each column
    for k in 1:9  # Each digit
        # Sum across columns (j) - row constraint
        @constraint(sudoku, sum(x[ind, j, k] for j in 1:9) == 1)
        # Sum across rows (i) - column constraint
        @constraint(sudoku, sum(x[i, ind, k] for i in 1:9) == 1)
    end
end
```
__Constraint 3__: we have the to enforce the constraint that each digit appears once in each of the nine ``3\times 3`` sub-grids. 
```julia
for i in 1:3:7
    for j in 1:3:7
        for k in 1:9
            # i is the top left row, j is the top left column.
            # We'll sum from i to i+2, e.g. i=4, r=4, 5, 6.
            @constraint(
                sudoku,
                sum(x[r, c, k] for r in i:(i+2), c in j:(j+2)) == 1
            )
        end
    end
end
```
"""

# ╔═╡ 420e8a1e-3756-41f2-987a-86acfc099456
begin
	radSolve = @bind radsolve Select(["Solve","Reset"],default="Reset")
	
	md"""
	$btnSkdo
	
	__Choose level of difficulty__
	$sdkoSelect

	__Action__
	$radSolve
	"""
end

# ╔═╡ 042a8b3b-fdc5-41f8-963e-3f51d147fe3b
# https://sugoku2.herokuapp.com/board?difficulty=easy
begin
	function getSudokuBoard(s)
		sd = Array{Int64, 2}(undef,9,9)
		for i in 1:9
			sd[i,:]=s[i]
		end
		sd
	end
	function printSudokuBoard(b)
		str=""
		for i in 1:9
			str*= (i % 3 == 0) ? "<tr style='border-bottom: 3px solid #ddd;'>" : "<tr>"
			for j in 1:9
				v = b[i,j]==0 ? "0" : "<b>$(b[i,j])</b>"
				str*= (j % 3 == 0) ? "<td style='border-right: 3px solid #ddd;'>" : "<td>"
				str*="$v</td>"
			end
			str*="</tr>"
		end
		"<table>$str</table>"
	end
	function printSudokuBoard(b,s)
		str=""
		for i in 1:9
			str*= (i % 3 == 0) ? "<tr style='border-bottom: 3px solid #ddd;'>" : "<tr>"
			for j in 1:9
				v = b[i,j]==0 ? "<span style='color:red;'>$(s[i,j])</span>" : "<b>$(b[i,j])</b>"
				str*= (j % 3 == 0) ? "<td style='border-right: 3px solid #ddd;'>" : "<td>"
				str*="$v</td>"
			end
			str*="</tr>"
		end
		"<table>$str</table>"
	end
	
end

# ╔═╡ 4a1396af-864f-4329-b833-6f2fea1d9309
function solveSudoku(initial_baord)
		sudoku = Model(HiGHS.Optimizer)
		set_silent(sudoku)
		@variable(sudoku,xx[i=1:9,j=1:9,k=1:9],Bin)
		
		for i in 1:9
		    for j in 1:9
		        # If the space isn't empty
		        if initial_baord[i, j] != 0
		            # Then the corresponding variable for that digit and location must
		            # be 1.
		            fix(xx[i, j, initial_baord[i, j]], 1; force = true)
		        end
		    end
		end
		for i in 1:9  # For each row
		    for j in 1:9  # and each column
		        # Sum across all the possible digits. One and only one of the digits
		        # can be in this cell, so the sum must be equal to one.
		        @constraint(sudoku, sum(xx[i, j, k] for k in 1:9) == 1)
		    end
		end
		for ind in 1:9  # Each row, OR each column
		    for k in 1:9  # Each digit
		        # Sum across columns (j) - row constraint
		        @constraint(sudoku, sum(xx[ind, j, k] for j in 1:9) == 1)
		        # Sum across rows (i) - column constraint
		        @constraint(sudoku, sum(xx[i, ind, k] for i in 1:9) == 1)
		    end
		end
		for i in 1:3:7
	    for j in 1:3:7
	        for k in 1:9
	            # i is the top left row, j is the top left column.
	            # We'll sum from i to i+2, e.g. i=4, r=4, 5, 6.
	            @constraint(
	                sudoku,
	                sum(xx[r, c, k] for r in i:(i+2), c in j:(j+2)) == 1
	            )
	        end
	    end
		end
		optimize!(sudoku)
		x_val = value.(xx);
		sol = zeros(Int, 9, 9)  # 9x9 matrix of integers
		for i in 1:9
		    for j in 1:9
		        for k in 1:9
		            # Integer programs are solved as a series of linear programs so the
		            # values might not be precisely 0 and 1. We can round them to
		            # the nearest integer to make it easier.
		            if round(Int, x_val[i, j, k]) == 1
		                sol[i, j] = k
		            end
		        end
		    end
		end
		sol
	end

# ╔═╡ 4260d8f7-2927-41ee-8c34-1e6ba37e4cbf
HTML("
<h2> References</h2>
<ol>
<li> <a name='wiki202208' href='https://en.wikipedia.org/wiki/Algebraic_modeling_language'>Algebraic modeling language</a>.(2022, August 28). In https://en.wikipedia.org/wiki/Algebraic_modeling_language</li>
<li>GAMS Documentation. <a name='gams' href='https://www.gams.com/latest/docs/'>https://www.gams.com/latest/docs/</a></li>
<li>Mosel Documentation. <a name='mosel' href='https://www.fico.com/fico-xpress-optimization/docs/latest/mosel/mosel_lang/dhtml/moselreflang.html'>https://www.fico.com/fico-xpress-optimization/docs/latest/mosel/mosel_lang/dhtml/moselreflang.html</a></li>
<li>LINGO Documentation. <a name='lingo' href='https://www.lindo.com/index.php/products/lingo-and-optimization-modeling'>https://www.lindo.com/index.php/products/lingo-and-optimization-modeling<a/></li>

<li>AIMMS Documentation. <a name='aimms' href='https://documentation.aimms.com/'>https://documentation.aimms.com/<a/></li>


<li>Pyomo Documentation. <a name='pyomo' href='http://www.pyomo.org/documentation'>http://www.pyomo.org/documentation<a/></li>


<li>JuMP Documentation. <a name='jump' href='https://jump.dev/JuMP.jl/v1.3.1/'>https://jump.dev/JuMP.jl/v1.3.1/<a/></li>

<li>
<a name='jump2017'> Dunning I, Huchette J, Lubin M (2017) JuMP: A Modeling Language for Mathematical Optimization. SIAM Review 59(2):295–320.</a>
</li>
<li><a name='cvx2006'>
M. Grant, S. Boyd, and Y. Ye, Disciplined convex programming, in Global Optimization:
From Theory to Implementation, Nonconvex Optimization and Its Application Series,
Springer, New York, 2006, pp. 155–210.</a>
</li>

<li><a name='cvx2014'>
M. Grant and S. Boyd, CVX: MATLAB Software for Disciplined Convex Programming,
Version 2.1, http://cvxr.com/cvx, 2014.</a>
</li>


<li><a name='fad2016'>
J. Revels, M. Lubin, and T. Papamarkou, Forward-Mode Automatic Differentiation in Julia, preprint, arXiv:1607.07892 [cs.MS], 2016; extended abstract presented at AD2016—7th International Conference on Algorithmic Differentiation, Oxford, UK.
</li>
<li>
<a name='ro2016'>I. Dunning, Advances in Robust and Adaptive Optimization: Algorithms, Software, and Insights, Ph.D. thesis, Massachusetts Institute of Technology, Cambridge, MA, 2016.</li>
</li>
</ol>
")

# ╔═╡ 5dcd42ad-19c1-4eee-9d89-a1bb5d6c67f2
begin
	struct SudokuBoard 
		board::Vector{Vector{Int64}}
	end
	function getLiveSudokuBoard(;live=:off,level="easy")
		
		if (live==:off)
		easy_boards =[
			[[0,0,0,0,0,8,0,0,5],[0,0,0,4,0,0,0,0,0],[0,6,0,1,0,0,2,0,7],[2,0,0,3,7,0,0,9,6],[3,5,6,0,9,1,4,7,0],[0,9,0,0,2,4,0,0,0],[0,4,5,7,1,0,0,0,0],[7,0,1,0,8,0,0,2,4],[9,8,2,5,0,3,0,6,0]],
			[[0,0,6,0,0,0,0,0,0],[0,2,0,5,0,8,0,0,9],[0,7,0,0,0,9,0,0,0],[0,0,0,0,0,0,0,0,0],[0,0,7,0,9,0,0,0,4],[0,0,9,4,2,7,3,0,0],[6,0,0,9,4,2,7,8,0],[7,0,2,6,8,5,9,1,0],[8,0,5,0,0,3,6,4,2]],
			[[0,5,7,1,0,0,8,0,0],[1,0,3,0,0,0,4,7,9],[0,0,0,2,3,7,0,0,6],[2,0,0,0,5,0,0,9,8],[0,6,0,0,9,0,2,0,0],[0,0,0,0,0,0,3,6,0],[0,0,0,6,2,4,9,8,7],[0,4,0,9,0,3,0,0,0],[9,0,0,8,1,0,6,3,4]]
		]
		medium_boards=[
			[[0,0,9,0,8,0,0,0,0],[1,2,3,4,0,7,0,0,0],[5,0,0,0,0,9,0,4,0],[0,1,4,3,0,0,0,9,8],[3,0,0,0,9,0,0,0,4],[0,9,0,0,0,0,0,0,0],[4,3,0,0,0,2,0,0,0],[6,7,0,9,0,0,0,0,3],[0,8,0,6,0,5,4,7,1]],
			[[0,5,0,0,9,0,0,0,0],[1,0,0,0,0,8,0,7,0],[7,8,9,0,0,0,0,0,0],[0,0,3,4,0,0,7,0,8],[4,6,5,0,0,0,0,2,3],[0,0,7,0,1,0,0,0,5],[5,0,0,0,4,0,0,0,0],[6,7,0,0,8,0,0,0,2],[0,3,0,5,2,7,6,0,1]],
			[[0,8,0,0,0,0,1,0,0],[0,0,0,0,5,0,0,8,0],[0,6,7,0,0,0,0,0,5],[0,0,0,0,0,0,7,9,8],[3,5,0,0,0,0,4,0,0],[0,0,0,0,0,0,0,0,0],[6,0,1,0,2,0,9,7,4],[8,0,2,0,7,0,0,0,0],[9,7,5,6,4,0,8,0,0]]
		]
		hard_boards=[
			[[0,0,0,0,6,0,2,0,0],[0,3,0,0,0,0,0,7,0],[0,0,0,2,0,0,0,0,0],[1,0,0,0,7,0,8,0,0],[4,5,0,0,0,0,0,1,0],[7,0,0,6,1,2,3,0,0],[3,0,0,5,2,6,9,0,0],[0,0,0,0,0,1,0,2,0],[0,0,0,8,0,0,0,6,1]],
			[[0,0,0,0,0,1,4,6,0],[0,0,0,0,0,6,0,8,0],[5,0,8,0,0,0,0,2,3],[2,1,3,5,0,0,6,0,0],[4,5,6,0,0,0,0,0,7],[0,0,0,0,2,0,0,0,0],[0,0,0,0,0,0,0,0,0],[6,0,0,0,0,4,8,5,0],[0,8,5,7,0,0,0,0,0]],
			[[7,0,0,0,0,2,6,3,5],[1,2,0,4,5,0,0,0,0],[0,0,0,0,0,8,1,0,4],[0,0,0,0,3,0,0,0,0],[3,4,6,0,8,0,2,0,0],[0,0,0,0,0,0,0,0,0],[0,3,1,0,0,7,0,0,8],[0,0,0,8,0,0,0,0,0],[0,7,0,0,0,0,0,1,0]]
		]
		boards = Dict(
			"easy" => easy_boards,
			"medium" => medium_boards,
			"hard" => hard_boards,
		)
			return rand(boards[level])
		else
		
		resp = HTTP.get("https://sugoku2.herokuapp.com/board?difficulty=$level")
		x3 = JSON3.read(resp.body, SudokuBoard)
		return x3.board
		end
	end
end

# ╔═╡ 1699ae56-2794-4277-b063-a69bb25ae840
begin
	data = [[3,0,4,0,0,0,0,5,0],[0,0,0,0,4,0,0,0,0],[0,0,9,1,5,7,0,3,4],[0,1,0,4,0,0,7,9,8],[4,0,6,0,0,9,0,0,3],[7,0,8,0,1,0,0,6,5],[0,3,0,8,0,0,9,4,0],[8,0,7,9,0,0,0,0,0],[0,0,0,0,0,0,3,8,0]]
	sdkoBtn
	initial_baord = getSudokuBoard(getLiveSudokuBoard(live=:off,level=lowercase(sdkoLevel)))
	
	HTML(printSudokuBoard(initial_baord))
end

# ╔═╡ 15fd06f8-1258-4edf-83c2-382d480f2d8f
begin
	if (radsolve=="Solve")
		sol = solveSudoku(initial_baord)
		HTML(printSudokuBoard(initial_baord,sol))
	else
		HTML(printSudokuBoard(initial_baord))
	end
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"
HiGHS = "87dc4568-4c63-4d18-b0c0-bb2238e4078b"
Ipopt = "b6b21f68-93f8-5de0-b562-5493be1d77c9"
JSON3 = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsPlots = "f3b207a7-027a-5e70-b257-86293d7955fd"

[compat]
HTTP = "~1.5.1"
HiGHS = "~1.2.0"
Ipopt = "~1.1.0"
JSON3 = "~1.11.1"
JuMP = "~1.4.0"
LaTeXStrings = "~1.3.0"
Plots = "~1.35.5"
PlutoUI = "~0.7.48"
StatsPlots = "~0.15.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "013dc2beba5a09216b8ceb023aaad6488cbe2693"

[[deps.ASL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6252039f98492252f9e47c312c8ffda0e3b9e78d"
uuid = "ae81ac8f-d209-56e5-92de-9978fef736f9"
version = "0.1.3+0"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "69f7020bd72f069c219b5e8c236c1fa90d2cb409"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.2.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "91ca22c4b8437da89b030f08d71db55a379ce958"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.3"

[[deps.Arpack_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "5ba6c757e8feccf03a1554dfaf3e26b3cfc7fd5e"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.1+1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "4c10eee4af024676200bc7752e536f858c6b8f93"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.1"

[[deps.BitFlags]]
git-tree-sha1 = "84259bb6172806304b9101094a7cc4bc6f56dbc6"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.5"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e7ff6cadf743c098e08fca25c91103ee4303c9bb"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.6"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "64df3da1d2a26f4de23871cd1b6482bb68092bd5"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.14.3"

[[deps.CodecBzip2]]
deps = ["Bzip2_jll", "Libdl", "TranscodingStreams"]
git-tree-sha1 = "2e62a725210ce3c3c2e1a3080190e7ca491f18d7"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.7.2"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "3ca828fe1b75fa84b021a7860bd039eaea84d2f2"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.3.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "46d2680e618f8abd007bce0c3026cb0c4a8f2032"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.12.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.DataValues]]
deps = ["DataValueInterfaces", "Dates"]
git-tree-sha1 = "d88a19299eba280a6d062e135a43f00323ae70bf"
uuid = "e7dc6d0d-1eca-5fa6-8ad6-5aecde8b7ea5"
version = "0.4.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "8b7a4d23e22f5d44883671da70865ca98f2ebf9d"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.12.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "04db820ebcfc1e053bd8cbb8d8bccf0ff3ead3f7"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.76"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "c36550cb29cbe373e95b3f40486b9a4148f89ffd"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.2"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "90630efff0894f8142308e334473eba54c433549"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.5.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "802bfc139833d2ba893dd9e62ba1767c88d708ae"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.5"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "187198a4ed8ccd7b5d99c41b69c679269ea2b2d4"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.32"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "00a9d4abadc05b9476e937a5557fcce476b9e547"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.69.5"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "bc9f7725571ddb4ab2c4bc74fa397c1c5ad08943"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.69.1+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "fb83fbe02fe57f2c068013aa94bcdf6760d3a7a7"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+1"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "a97d47758e933cd5fe5ea181d178936a9fc60427"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.5.1"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HiGHS]]
deps = ["HiGHS_jll", "MathOptInterface", "SparseArrays"]
git-tree-sha1 = "d40a9e8db6438481915261a378fc2c8ca70bb63a"
uuid = "87dc4568-4c63-4d18-b0c0-bb2238e4078b"
version = "1.2.0"

[[deps.HiGHS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3e24275666fcc2e24d1a58a9f02acd9d2e23d3a"
uuid = "8fd58aa0-07eb-5a78-9b36-339c94fd15ea"
version = "1.3.0+0"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions", "Test"]
git-tree-sha1 = "709d864e3ed6e3545230601f94e11ebc65994641"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.11"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "842dd89a6cb75e02e85fdd75c760cdc43f5d6863"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.6"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.Ipopt]]
deps = ["Ipopt_jll", "MathOptInterface"]
git-tree-sha1 = "14a305ededd75330246aaa0380130561d8924120"
uuid = "b6b21f68-93f8-5de0-b562-5493be1d77c9"
version = "1.1.0"

[[deps.Ipopt_jll]]
deps = ["ASL_jll", "Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "MUMPS_seq_jll", "OpenBLAS32_jll", "Pkg"]
git-tree-sha1 = "e3e202237d93f18856b6ff1016166b0f172a49a8"
uuid = "9cc047cb-c261-5740-88fc-0cf96f7bdcc7"
version = "300.1400.400+0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JSON3]]
deps = ["Dates", "Mmap", "Parsers", "SnoopPrecompile", "StructTypes", "UUIDs"]
git-tree-sha1 = "65edf3850efb9cb4ca3b0bf488e29c6c38a23d2d"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.11.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.JuMP]]
deps = ["LinearAlgebra", "MathOptInterface", "MutableArithmetics", "OrderedCollections", "Printf", "SparseArrays"]
git-tree-sha1 = "9a57156b97ed7821493c9c0a65f5b72710b38cf7"
uuid = "4076af6c-e467-56ae-b986-b466b2749572"
version = "1.4.0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "9816b296736292a80b9a3200eb7fbb57aaa3917a"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.5"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "ab9aa169d2160129beb241cb2750ca499b4e90e9"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.17"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "94d9c52ca447e23eac0c0f074effbcd38830deb5"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.18"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "5d4d2d9904227b8bd66386c1138cf4d5ffa826bf"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "0.4.9"

[[deps.METIS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "1d31872bb9c5e7ec1f618e8c4a56c8b0d9bddc7e"
uuid = "d00139f3-1899-568f-a2f0-47f597d42d70"
version = "5.1.1+0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "2ce8695e1e699b68702c03402672a69f54b8aca9"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.2.0+0"

[[deps.MUMPS_seq_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "METIS_jll", "OpenBLAS32_jll", "Pkg"]
git-tree-sha1 = "29de2841fa5aefe615dea179fcde48bb87b58f57"
uuid = "d7ed1dd3-d0ae-5e8e-bfb4-87a502085b8d"
version = "5.4.1+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MathOptInterface]]
deps = ["BenchmarkTools", "CodecBzip2", "CodecZlib", "DataStructures", "ForwardDiff", "JSON", "LinearAlgebra", "MutableArithmetics", "NaNMath", "OrderedCollections", "Printf", "SparseArrays", "SpecialFunctions", "Test", "Unicode"]
git-tree-sha1 = "2284cb18c8670fd5c57ad010ce9bd4e2901692d2"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "1.8.2"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.MultivariateStats]]
deps = ["Arpack", "LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI", "StatsBase"]
git-tree-sha1 = "efe9c8ecab7a6311d4b91568bd6c88897822fabe"
uuid = "6f286f6a-111f-5878-ab1e-185364afe411"
version = "0.10.0"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "1d57a7dc42d563ad6b5e95d7a8aebd550e5162c0"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.0.5"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "440165bf08bc500b8fe4a7be2dc83271a00c0716"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "5a9ea4b9430d511980c01e9f7173739595bbd335"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.2"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "f71d8950b724e9ff6110fc948dff5a329f901d64"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS32_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c6c2ed4b7acd2137b878eb96c68e63b76199d0f"
uuid = "656ef2d0-ae68-5445-9ca0-591084a874a2"
version = "0.3.17+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "3c3c4a401d267b04942545b1e964a20279587fd7"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.3.0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e60321e3f2616584ff98f0a4f18d98ae6f89bbb3"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.17+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "cf494dca75a69712a72b80bc48f59dcf3dea63ec"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.16"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "6c01a9b494f6d2a9fc180a08b182fcb06f0958a0"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.4.2"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "21303256d239f6b484977314674aef4bb1fe4420"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SnoopPrecompile", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "0a56829d264eb1bc910cf7c39ac008b5bcb5a0d9"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.35.5"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "efc140104e6d0ae3e7e30d56c98c4a927154d684"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.48"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "c6c0f690d0cc7caddb74cef7aa847b824a16b256"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "97aa253e65b784fd13e83774cadc95b38011d734"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.6.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "d12e612bba40d189cead6ff857ddb67bd2e6a387"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase", "SnoopPrecompile"]
git-tree-sha1 = "9b1c0c8e9188950e66fc28f40bfe0f8aac311fe0"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.7"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "efd23b378ea5f2db53a55ae53d3133de4e080aa9"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.16"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SnoopPrecompile]]
git-tree-sha1 = "f604441450a3c0569830946e5b33b78c928e1a85"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "f86b3a049e5d05227b10e15dbb315c5b90f14988"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.9"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "5783b877201a82fc0014cbf381e7e6eb130473a4"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.0.1"

[[deps.StatsPlots]]
deps = ["AbstractFFTs", "Clustering", "DataStructures", "DataValues", "Distributions", "Interpolations", "KernelDensity", "LinearAlgebra", "MultivariateStats", "NaNMath", "Observables", "Plots", "RecipesBase", "RecipesPipeline", "Reexport", "StatsBase", "TableOperations", "Tables", "Widgets"]
git-tree-sha1 = "e0d5bc26226ab1b7648278169858adcfbd861780"
uuid = "f3b207a7-027a-5e70-b257-86293d7955fd"
version = "0.15.4"

[[deps.StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "ca4bccb03acf9faaf4137a9abc1881ed1841aa70"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.10.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableOperations]]
deps = ["SentinelArrays", "Tables", "Test"]
git-tree-sha1 = "e383c87cf2a1dc41fa30c093b2a19877c83e1bc1"
uuid = "ab02a1b2-a7df-11e8-156e-fb1833f50b87"
version = "1.2.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "8a75929dcd3c38611db2f8d08546decb514fcadf"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.9"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "e59ecc5a41b000fa94423a578d29290c7266fc10"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "fcdae142c1cfc7d89de2d11e08721d0f2f86c98a"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.6"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╟─c154e809-bee8-444f-921f-32100509519e
# ╟─07d5b126-4150-11ed-24e1-eb98b921c8c6
# ╟─d716ec58-253d-4267-9497-7497fe639c6a
# ╟─bcfd72c1-c3f5-4cc6-96dd-b9d292d78c7f
# ╟─f111a4cf-7f17-429b-9bb0-290fe484e2d0
# ╟─b193fb26-5736-4b37-9c8c-6a44072b4a6d
# ╟─01488a29-905c-4ea4-8fee-d0c49b81f1e0
# ╟─3327bf8e-d5d3-44f5-85c1-c1bde64abcc3
# ╟─111e9ba3-7a25-4a27-90f3-fe358d1a66a5
# ╟─5b57a081-d9e7-42d3-b579-8cf579d851bf
# ╟─0e4cdd43-3fc0-467f-b004-250fbaf53ac0
# ╟─85c4eb73-c2b2-40e3-85f8-09d5639f5d37
# ╟─b51da029-f797-45d9-9309-c18050474ca6
# ╟─ab16ed26-4aa7-40d7-8013-a2ed51352624
# ╟─ab450756-5cfc-4df3-9bb2-ad0b9da8c948
# ╟─46e87d72-d899-43dd-9c81-7eb01200e6e7
# ╟─c57bcac2-9cb7-4c6d-ac7f-3436d3cbcc6f
# ╟─43a9716e-d735-47d1-b209-77754bad82ab
# ╟─78e26d05-cd83-4029-a300-b21a38f7d724
# ╟─738f6aa6-7780-41c4-89ef-948946ba81b6
# ╟─197a354d-e466-4aeb-8590-b176f6fc96bd
# ╟─38573590-c70a-413c-a141-8733b0ee1a6e
# ╟─53357188-bb77-41aa-b7a5-3dfe7ff4739e
# ╟─57855292-adf7-4654-9583-c76d241f055b
# ╟─745f32fb-1f54-4cd2-b6ff-73df903165a3
# ╟─c769c0ab-9f28-4939-812e-4329a8242d14
# ╟─76a3d98d-c700-476e-92a7-a4d5cbf686cf
# ╟─c0c3a7dd-5a54-4966-b980-c786a6861ede
# ╟─d6fa6332-d7c6-4283-906e-b25b66fa6067
# ╟─296d6c22-0ac1-4003-b93e-6d0261402e3c
# ╟─1699ae56-2794-4277-b063-a69bb25ae840
# ╟─06c991d3-c550-44b9-99b1-563c81088066
# ╟─420e8a1e-3756-41f2-987a-86acfc099456
# ╠═15fd06f8-1258-4edf-83c2-382d480f2d8f
# ╟─042a8b3b-fdc5-41f8-963e-3f51d147fe3b
# ╟─4a1396af-864f-4329-b833-6f2fea1d9309
# ╟─4260d8f7-2927-41ee-8c34-1e6ba37e4cbf
# ╟─5dcd42ad-19c1-4eee-9d89-a1bb5d6c67f2
# ╠═a563c4b9-01d0-4411-97ae-8019414330b1
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
