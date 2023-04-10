md"""

# First steps

## Comments

Comments are non-coding part in the source code. Although the compiler does not read comments, comments are important notes for humans, making the code more readable (hopefully).

```julia
# One line comments

#=
Multiline comments
=#
```

## Variables and expressions

The assignment operator `=` binds a name to a piece of data.

To see the data content:

- `@show expression` to show the expression and its result
- [`println(data)`](https://docs.julialang.org/en/v1/base/io-network/#Base.println): good old print function.
- [`@printf`](https://docs.julialang.org/en/v1/stdlib/Printf/): C-like formatted output.
- Inline display: content of the last expression. A semicolon `;` at the end will suppress inline output.

"""

# An integer (64 bit)
x = 1

# A (64-bit) floating-point number
y = 1.0

# Also a (64-bit) floating-point
y = 1.

# A 32-bit floating-point number. less accurate but calculates faster. Often used in GPU computing.
y = 1.0f0

# Complex number
z = 3 + 5im

# Unicode names: "\alpha<tab>"
α = 3.74

# Strings (Text) are surrounded by double quotes. NOT SINGLE QUOTES!
s = "Julia"

# Characters are surrounded by single quotes
c = ' '

# Fractions (rational numbers)
station = 9 + 3//4

# Constants will emit a warning if you try to change it after its creation
const theUltimateAnswer = 42

#---
println("Hello World")

#---
println("Hello ", s)

#---
println(1, 2, 3)

# @show will print x = val
@show x	1+1 2-2 3*3;

#---
@show typeof(x) typeof(y) typeof(c)  typeof(s) typeof(z) typeof(1//2);

# `convert(T,x)`` converts x to type T
typeof(convert(Float64, x))

# There is also `Type(x)`
typeof(Float64(x))

# Or this convenience function
typeof(float(x))

#===
## Compound expressions

- A [begin block](https://docs.julialang.org/en/v1/base/base/#begin) `begin` ... `end` squashes multiple expressions into one.
- A [let block](https://docs.julialang.org/en/v1/manual/variables-and-scoping/#Let-Blocks) `let` ... `end` is similar to a begin block but variables inside will be discarded outside.
===#

# a1 and a2 are available after begin block ends
begin
	a1 = 2
	a2 = 35
	a1 = a1 * a2
end

# x, y, z are NOT available after let block ends
let
	x = 1
	y = 2
	z = 3

	x + y * z
end

#===
## Math expressions

Julia supports basic arithmatic operations and essential math functions by default.

### Basic arithmetic
===#

# Multiple assignment
a, b = 2, 3

# Addition
a + b

# Subtraction
a - b

# Multiplication
a * b

# Division
b / a

# Fraction
b // a

# integer division: `\div<tab>`, the same as `div(a, b)`
a ÷ b

# Modulus
b % a

# Power
a^b

#===
## Comparison

Returns a Boolean value (`true` / `false`)
===#

a, b = 2, 3

#---
a < 1

#---
b > 2

#---
a <= b

#---
a != b

#---
a == b + 1

# Chained comparisons are supported
1 < a <= b

# Approximation operator `\approx<TAB>` for floating point number equivalence
1e10 + 1.0 ≈ 1e10

# The same as
isapprox(1e10 + 1.0, 1e10)

# ## Math functions

# `\pi<TAB>`
sin(0.5*π)

# More precise
sinpi(1//2)

#---
cos(0.5*π)

# More precise
cospi(1//2)

# `\sqrt<TAB>`
sqrt(π) ≈ √π

# Natual log
log(10)

# Common log
log10(10)

# Natural exponant
exp(-5)

# expm1(x) is more accurate that exp(x) - 1 when x is very close to zero
exp(1e-12) - 1, expm1(1e-12)

#===

## Strings

- [Think Julia ch8](https://benlauwens.github.io/ThinkJulia.jl/latest/book.html#chap08)
- [Introduction_to_Julia_tutorials](https://nbviewer.jupyter.org/github/xorJane/Introduction_to_Julia_tutorials/blob/master/02.Strings.ipynb)

A `string` is a sequence of characters.

- `" ... "` for one line strings.
- Three double quotes surround multiline string.
- `str1*str2*...` to concatenate strings
- `string(str1, str2, ...)` to convert the data (if neede) and make a string.
- `^` to repeat a string: `str^3` to repeat `str` three times.
- `[idx]` to access an individual character.
- `$` to insert (or interpolate) a value into a string.


Although `string(x, y)` looks less streamlined, it is generally faster than interpolation `$` or concatenation `*` and is most recommended.

===#

"I am a string."

#---

"""
I am a multiline
string.

Hello Julia!
"""

# A character is different from a string
'a' == "a"

# How to insert contents into a string
str1 = "BEBI"
str2 = "5009"
string("The class is ", str1, '-', str2, ".")

# Use string interpolation
"The class is $str1-$str2."

# concat string using `*`
str1*"-"*str2

#===

## Control flow

- [control flow](https://docs.julialang.org/en/v1/manual/control-flow/)
- [functions](https://docs.julialang.org/en/v1/manual/functions/)


Julia programs are able to run nonlinearly by controlling its execution flow.

### Conditional statements

- The `elseif` and `else` blocks are optional. `if` and `end` are mandatory.
- `if` blocks return a value. Capturing the value is optional.
- `if` blocks are "leaky", i.e. they do not introduce a local scope. (The same as Python)
- Only boolean (true / false) could be evaluated in `if` blocks. Using other types (e.g. Int) will generate an error.

#### ifelse function

All the arguments are evaluated first in `ifelse(cond, tvalue, fvalue)`.


#### short cicuit evaluation

`&&` (logical and) and `||` (logical or) operators support [short circuit evaluation](
Short-circuit evaluation - Wikipedia
https://en.wikipedia.org › wiki › Short-circuit_evaluation).

- In the expression `a && b`, the subexpression `b` is only evaluated if `a` evaluates to true.
- In the expression `a || b`, the subexpression `b` is only evaluated if `a` evaluates to false.

===#

# && evaluates and returns the second argument if the first is true
true && println("Hi")

# && otherwise returns false
false && println("Hi")

# `if` block has return value(s)
let score = 10
    response = if 80 < score <= 100
        "Good"
    elseif 60 < score <= 80
        "Okay"
    else
        "Uh-Oh"
    end
    response
end

# Ternary operator
1 > 2 ? "True" : "False"

#===
### Loops

Loops are repeated evaluations of a code block.

- `while` loops are often related to a predicate.
- `for` loops are often related to a sequence.


- `break`: exit the loop immediately.
- `continue`: move on to the next item / evaluation immediately.
===#

# Hailstone sequence (3n+1 problem) in a while loop
let n = 1025489
    step = 0
    while n > 1
        if iseven(n)
            n ÷= 2
        else
            n = 3n + 1
        end
        step += 1
    end
    step
end

# Summation
let n = 100
    s = 0
    for i in 1:n
        s += i
    end
    s
end

# For loop
for x in 1:9
    if x == 5
        continue ## jump to line #2
    elseif x >=8
        break    ## jump to line #9
    end
    println(x, "^2 = ", x^2)
end

# Use enumerate(seq) to get a pair of idx and element
for (i, x) in enumerate(10:-1:1)
    println("xs[", i, "] = ", x)
end

# Multiple nested for loops can be combined into a single outer loop, forming the cartesian product of its iterables.
for i = 'x':'z', j = '1':'3'
    println(i, j)
end

#===
## Functions

> In Julia, a function is an object that maps a tuple of argument values to a return value.
> [Julia docs](https://docs.julialang.org/en/v1/manual/functions/)


Functions facilitate:
- Code reuse and encapsulation.
- Specializations ([Methods](https://docs.julialang.org/en/v1/manual/methods/))

- Functions are first-class objects and can be passed into a higher-order function.
- The arguments are "passed-by-sharing". Modifications to mutable argument values (such as `Arrays`) will be reflected to the caller. (Similar to Python)
- Functions that will update the arguments are named with a bang `!` by convention. (e.g. sort(arr) vs sort!(arr))
- Often only the scalar version of a function is required; for vector element-wise operations, there are broadcast (dot) syntax.
- You can write multiple function with the same name provided they have different parameter lists. Julia will choose the most apporpriate one for you.

### Standard form
===#

"Mechaelis-Menton function"  ## Function documentations
function mm(x, k)            ## function name and parameter list
  result = x / (x +k)        ## Doing calculations
  return result              ## return statement is optional
end

#---
mm(1, 0.5)

# ### One-liner form

f(x, y) = x + y
f(1, 2)

# And you can also reuse previously-defined functions
mm(x) = mm(x, one(x))
mm(1)

#===
### Anonymous functions
Often used with high-order functions e.g. `map()`
===#

g = x->x^2
g(3)

#---
map(x->x^2, 1:3)

# Use [`do` block](https://docs.julialang.org/en/v1/manual/functions/#Do-Block-Syntax-for-Function-Arguments) for long anonymous functions.

val = rand(-6:6, 10)

map(val) do x
    res = if x < 0 && iseven(x)
        zero(x)
    elseif x == 0
        one(x)
    else
        x
    end
    res
end

#===

### Optional arguments

[Optional (positional) arguments](https://docs.julialang.org/en/v1/manual/functions/#Optional-Arguments) are listed after mandatory ones.

```julia
function func(a, b, c=1)

end
```

And they are called with `func(a, b)` or `func(a, b, 3)`

### Keyword arguments

[Keyword arguments](https://docs.julialang.org/en/v1/manual/functions/#Keyword-Arguments) are listed after `;`

```julia
function plot(x, y; style="solid", width=1, color="black")
    ...
end
```

And they are called with `plot(x, y, width=2)` or `plot(x, y; width=2)`

===#

args_kwargs(args...; kwargs...) = (args, kwargs)  ## mind the semicolon ;

args_kwargs(1, 2, 3; a=4, b=5.0, c="Hello")

#===

### See also

- [Compositing functions and pipes](https://docs.julialang.org/en/v1/manual/functions/#Function-composition-and-piping)
- [Variable argument (vararg) functions](https://docs.julialang.org/en/v1/manual/functions/#Varargs-Functions)
- [Argument destructuring](https://docs.julialang.org/en/v1/manual/functions/#Argument-destructuring)
- [Scope of variables](https://docs.julialang.org/en/v1/manual/variables-and-scoping/)

===#

#===

## Collections, broadcasting, and Methods
Using built-in collections is the simplest way to group and organize data.

The values in a `immutable` collection cannot be updated after its creation, while in a `mutable` collection can.

The elements in `sequential` collections are accessed by integer indices, while those in `associative` collection are accessed by keys.

### Sequential collections

General rules for sequential collections:

- 1-based indexing, as in R, MATLAB, and Fortran.
- Elements are accessed by an integer index `seq[i]` or an integer range `seq[1:2:end-1]`.
- `length(seq)` returns the total size
- Splat operator `...` passes the inner contents in the collection as positional function arguments.
- Dot syntax (e.g. a .+ b) performs element-wise  / broadcasting operations.

#### Ranges

`start[:step]:end`

- Immuatable
- Sequencial
- Eheap
- Evenly-spaced numerical sequences

===#

# A simple range
1:2:10

# Length of a sequence
length(1:2:10)

# Show its content
dump(1:2:10)

# Explicit range function
range(1, 10, length=10)

# linspace() equivalent
LinRange(1, 10, 10)

# Pick elements by a range of indices
(1:10)[3:end]

#===

#### Tuples

- immutable
- sequential collections
- efficient for heterogenous data of various types
- stack-allocated

===#

tuple(1, 'a', 3.14)

# Usually written as
(1, 'a', 3.14)

# Accessing elements
t1 = (1, 2, 3)
t1[1]

#---
t2 = (1, 'a', 3.14)
dump(t2)

# Merging multiple tuples using the splat (...) operator
tuple(t1..., t2...)

## Tuples could be used to swap elements
let x = 1, y = 2, z = 3
	x, y, z = z, x, y
	@show x, y, z
end;

# Tuple can return multiple values from a function

neighbors(x) = x+1, x-1
neighbors(0)

#---
extrema([1, 5, 6, 7, -1, -3, 0])

#---
sincospi(1//2)

#===

#### Arrays

[Arrays in Julia docs](https://docs.julialang.org/en/v1/manual/arrays/)

`[seq...]` / `collect(seq)`

Arryas are the bread and butter for scientific computing, similar to numpy's ndarrays.

- Column-major (Fortran style) rather than row-major (C and numpy style)
- Assignments and updating may cuase unwanted editing due to memory sharing.

Some useful functions for arrays:

- `length(A)` the number of elements in A
- `ndims(A)` the number of dimensions of A
- `size(A)` a tuple containing the dimensions of A
- `size(A,n)` the size of A along dimension n
- `eachindex(A)` an efficient iterator for visiting each position in A

===#

# 1D array (column vector)
x = [5, 6, 7]

# np.arange() equivalent
collect(1:10)

# Array with all zeroes
zeros(2, 5, 2)

# Array with all ones
ones(2, 5)

# Uninitialized array with the same data type and dims as x
similar(x)

# np.zeros_like()
zero(x)

# Array of random numbers
rand(1:6, 10)

#---
rand(1:6, 2, 2)

# Reshape an array
reshape(1:12, 3, 4)

# Reshape A to an (1D) vector
vec(rand(1:6, 2, 2))

# repeat the array 3x2 times
repeat([1 2; 3 4], 3, 2)

# comprehension for 1D array
[i^2 for i in 1:10 if i >= 5]

# 2D comprehension for a 2x3 array
[x * y for x in 1:2, y in 1:3]

# casting comprehension result element type to Float64
Float64[x^2 for x in 1:4]

# This is a 1-element tuple containing a vector
tuple([1,2,3])

# How to convert vector to tuple
Tuple([1,2,3])

# 2D array (matrix)
# space is a shorthand for hcat()
# semicolone is a shorthand for vcat()
A = [1 2 3;
     4 5 6]

# Accessing elements
A[1, 2]

# Accessing a range of elements
A[1:2, 2:3]

# Array total length
length(A)

#---
axes(A)

#---
size(A)

#---
ndims(A)

#---
transpose(A)

# (Conjugate transpose) Adjoint
A'

# Matrix-vector multiplication
b = A * x

# Find x for Ax = b, using left division operator `/`
A\b ≈ x

# Flatten A to an (1D) vector
vec(A)

# Arrays are mutable (i.e. you can update the contents) objects
# You should make a copy if you want the original one intact
A[1, 1] = 0
A

#===

### Associative collections

- `d[key]` accesses values by keys
- `d[key] = value` sets a key-value pair for a mutable dictionary.
- `delete!(d, key)` deletes the kay (and its partner) from a mutable dictionary.
- `keys(d)` returns a series of keys
- `values(d)` returns a series of values
- `pairs(d)` returns a series of (key => value) pairs
- `merge(d1, d2, ...)` return combinations of several dicts. `merge!(d1, d2, ...)` combine several dicts and update the first one.
- `get(d, key, default)` returns the value stored for the given key, or the given default value if no mapping for the key is present.

#### Named tuples

`Namedtuple`s are tuples with key-value pairs.

===#

nt = (a=1, b=2, c=4)

#---
nt[1]

#---
nt.a == nt[:a] == nt[1]

# How to fill a named tuple elegantly

a = 1
b = 2
c = 3

nt = (; a, b, c)

#===
#### Dictionaries

[Dictionaries](https://docs.julialang.org/en/v1/base/collections/#Dictionaries) are mutable mappings of `key => value`.
===#

eng2sp = Dict("one" => "uno", "two" => "dos", "three" => "tres")

#---
eng2sp["two"]

#---
eng2sp["five"] = "cinco"

#---
keys(eng2sp)

#---
values(eng2sp)

#---
get(eng2sp, "one", "N/A")

#---
get(eng2sp, "four", "N/A")

#---
haskey(eng2sp, "one")

# Elements are not ordered
for (k ,v) in eng2sp
    println(k, " => ", v)
end

# Creating a dict from an arrya of tuples
Dict([('a', 1), ('c', 3), ('b', 2)])

# Creating a Dict via a generator (similar to comprehensions)
Dict(i => i^2 for i = 1:10)

#---
Dict(zip("abc", 1:3))

#===

## Broadcasting (Dot) syntax

[Broadcasting](https://docs.julialang.org/en/v1/manual/arrays/#Broadcasting) turns scalar operations into vector ones.

===#

[1, 2, 3] .+ [4, 5, 6]

#---
[1, 2, 3] .+ 4

# Element-wise operation
sinpi.([0.5, 1.0, 1.5, 2.0])

# Create a dictionary with a list of keys and a list of values
ks = (:a, :b, :c)
vs = (1, 2, 3)
Dict(ks .=> vs)

# How to do logspace() in Julia
exp10.(LinRange(-3.0, 3.0, 50))

# Make a 9*9 multiplication table
collect(1:9) .* transpose(collect(1:9))

#===
## Custom data structures

https://docs.julialang.org/en/v1/manual/types/#Composite-Types
===#

# struct or mutable struct
struct Point
    x
    y
end

# Define a default constructor
Point() = Point(0.0, 0.0)

#---
p1 = Point(1.0, 2.0)
p2 = Point(-3.0, 2.0)

# Define a method for our cutsom type
add(a::Point, b::Point) = Point(a.x + b.x, a.y + b.y)

#---
add(p1, p2)

#===
## Methods

[Methods in Julia docs](https://docs.julialang.org/en/v1/manual/methods/#Methods)

You can overload the same function with different argument types/numbers. Julia will try to find the right function for the argument type(s).
===#

func(a::Int) = a + 2
func(a::AbstractFloat) = a/2
func(a::Rational) = a/11
func(a::Complex) = sqrt(a)
func(a, b::String) = "$a, $b"

#---
func(1)

#---
func(3.0)

#---
func(33//4)

#---
func(-2 + 0im)

#---
func(true, "it just works and compiles down to optimal code")

# ## Runtime information
import Pkg
Pkg.status()

#---
import InteractiveUtils
InteractiveUtils.versioninfo()
