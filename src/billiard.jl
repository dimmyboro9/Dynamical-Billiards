"""

Draw and simulate motion of a particle confined to a billiard with boundary
formed by line segments and circular arcs.

"""
module Billiard

using LinearAlgebra
using Luxor
using CairoMakie

export Arc, Segment, Table
export simulate!, draw_path, draw_poincare

const K_BORDER = 5000

#
# Types
#

"""

    Ray(pos, velocity)

A ray starting at the point "pos" and directed in the direction of the "velocity" vector.

During the search for the intersection point, I consider the ray as a line, expressing it by the equation y = kx + b or x = b (for lines parallel to the y-axis). 
In the first case, `k` is the slope of the line, the tangent of the angle between the line and the positive x-axis, 
and `b` is the y-intercept, the point where the line intersects the y-axis. 
In the second case, `k` is undefined, and `b` is a specific number determining the constant value of x.

"""
struct Ray
    pos::Vector{Float64}
    velocity::Vector{Float64}
    k::Float64
    b::Float64

    function Ray(p::Vector{Float64}, v::Vector{Float64})
        k = isapprox(v[1], 0) ? NaN : v[2] / v[1]
        b = isnan(k) ? p[1] : p[2] - k * p[1]
        new(p, v, k, b)
    end
end

abstract type Obstacle end

"""

    Segment(A, B)

Obstacle formed by a line segment connecting points `A` and `B`. 

During the search for the intersection point, I consider the segment as a line, expressing it by the equation ``y = kx + b`` or ``x = b`` (for lines parallel to the y-axis). 
In the first case, `k` is the slope of the line, the tangent of the angle between the line and the positive x-axis, 
and `b` is the y-intercept, the point where the line intersects the y-axis. 
In the second case, `k` is undefined, and b is a specific number determining the constant value of x.

The vector `s` represents the directional vector, indicating the direction of the line on which the line segment lies. 
The normal vector `n` is any vector perpendicular to `s` with a length equals 1.

"""
struct Segment <: Obstacle
    A::Vector{Float64}
    B::Vector{Float64}
    s::Vector{Float64}
    n::Vector{Float64}
    k::Float64
    b::Float64

    function Segment(A::Vector{Float64}, B::Vector{Float64})
        (length(A) != 2 || length(B) != 2 || isapprox(A, B)) && throw(ArgumentError("Invalid Segment arguments. Please provide two distinct vectors, each with exactly two elements."))
        
        k = isapprox(A[1], B[1]) ? NaN : (A[2] - B[2]) / (A[1] - B[1])
        b = isnan(k) ? A[1] : A[2] - k * A[1]
        
        s = (A - B) / norm(A - B)
        n = isapprox(s[2], 0.0) ? [s[2], -s[1]] : [-s[2], s[1]]
        
        new(A, B, s, n, k, b)
    end
end

@doc raw"""

    Arc(C, R, α, β)

Obstacle formed by a circular arc defined by the circle of radius `R` and
center `C`, and specified by the polar coordinates `α` and `β` (measured counter clock-wise with respect to `R` and with the angle running from ``-π`` to ``π``, the direction of the positive ``x`` axis corresponds to
the angle ``0``).

More precisely, the circular arc is parametrized by ``x = R \cos(t)`` and
``y = R \sin(t)` for ``t ∈ ⟨α, β⟩`` where ``α < β``.

"""
struct Arc <: Obstacle
    C::Vector{Float64}
    R::Float64
    α::Float64
    β::Float64

    function Arc(C::Vector{Float64}, R::Float64, α::Union{Float64, Irrational{:π}}, β::Union{Float64, Irrational{:π}})
        # In my code, everywhere >= is defined with > and isapprox(). All this is done because we work with machine numbers.
        # a >= b := a > b || isapprox(a, b);
        # a <= b := a < b || isapprox(a, b).  
        (length(C) != 2 || R < 0.0 || β < α || isapprox(β, α)) && throw(ArgumentError("Invalid Arc arguments. Please provide the center of the arc as a vector of two elements, a non-negative radius, and two angles α and β, where α is less than β."))
        α = mod2pi(Float64(α))
        β = mod2pi(Float64(β))
        
        α = α > π || isapprox(α, π) ? α - 2 * π : α
        β = β < π || isapprox(β, π) ? β : β - 2 * π
        
        new(C, R, α, β)
    end
end

"""

    Table(obstacles)

The billiard table is formed by number of obstacles.

"""
struct Table
    obstacles::Vector{<:Obstacle}

    function Table(obstacles::Vector{<:Obstacle})
        new(obstacles)
    end
end

#
# Functions and methods
#

"""

    satisfies_constraints(obstacle/ray, point)

A bunch of functions that check whether the found point lies on a ray, line or arc
"""
function satisfies_constraints(ray::Ray, point::Vector{Float64})
    point_ray = point - ray.pos
    
    if isapprox(ray.velocity[1], 0.0) 
        # Here I have a question. I don't know why but in Julia there is a difference between:
        #
        # isapprox(82.81000313634038, 82.8100031363404) -> returns true
        # isapprox(82.81000313634038 - 82.8100031363404, 0.0) -> returns false
        # isapprox(82.81000313634038 / 82.8100031363404, 1.0) -> returns true
        #
        # I would be grateful if you could explain to me why it works this way
        #
        # I'm guessing that the point is that there is a more machine numbers around zero, but I'm not sure about it
        
        (isapprox(point[2], ray.pos[2]) || point_ray[2] / ray.velocity[2] < 0.0) && return false
    else
        (isapprox(point[1], ray.pos[1]) || point_ray[1] / ray.velocity[1] < 0.0) && return false
    end
    
    return true
end

function satisfies_constraints(segment::Segment, point::Vector{Float64})
    point_a = point - segment.A
    b_a = segment.B - segment.A
    
    if isapprox(b_a[1], 0.0)
        t = point_a[2] / b_a[2]
        ((! isapprox(point[2], segment.A[2]) && t < 0.0) || (! isapprox(t, 1.0) && t > 1.0)) && return false
    else
        t = point_a[1] / b_a[1]
        ((! isapprox(point[1], segment.A[1]) && t < 0.0) || (! isapprox(t, 1.0) && t > 1.0)) && return false
    end
    
    return true
end

function satisfies_constraints(arc::Arc, point::Vector{Float64})
    point_relative_C = point - arc.C
    point_angle = (point_relative_C[2] > 0.0 ? 1.0 : -1.0) * acos(dot(point_relative_C, [1, 0]) / norm(point_relative_C))
    
    on_arc = arc.α < arc.β ? 
            (point_angle > arc.α || isapprox(point_angle, arc.α)) && (point_angle < arc.β || isapprox(point_angle, arc.β)) : 
            (point_angle > arc.α || isapprox(point_angle, arc.α)) || (point_angle < arc.β || isapprox(point_angle, arc.β)) 
    
    return on_arc 
end

"""

    solve_quadratic_equation(a, b, c)

Functions that implement the search for the roots of a quadratic equation
"""
function solve_quadratic_equation(a, b, c)        
    discriminant = b^2 - 4 * a * c
    
    if isapprox(discriminant, 0.0)
        return 1, [-b / (2*a)]
    elseif discriminant > 0.0
        discr_root = √discriminant
        return 2, [(-b + discr_root) / (2*a), (-b - discr_root) / (2*a)]
    else 
        return 0, NaN
    end
end

"""

    transition_qe(k, b, m, n, R)

A transition function to find the intersection of a line with a circle, i.e. to solve a quadratic equation. 
If the coefficient k is too large, it switches to the BigFloat type in order not to lose the accuracy of calculations due to working with machine numbers.
"""
function transition_qe(k, b, m, n, R)
    if abs(k) > K_BORDER
        k, b = BigFloat(k), BigFloat(b)
        num_of_roots, roots = solve_quadratic_equation(k^2 + 1, 2 * (k * (b - n) - m), m^2 + (b - n)^2 - R^2)
        return num_of_roots, Float64.(roots)
    end
    return solve_quadratic_equation(k^2 + 1, 2 * (k * (b - n) - m), m^2 + (b - n)^2 - R^2)
end

"""

    find_intersection(obstacle, obstacle/ray)

A bunch of functions that looks for an intersection point between different obstacles and ray
"""
function find_intersection(segment::Segment, ray::Ray)
    (isapprox(segment.k, ray.k) || (isnan(segment.k) && isnan(ray.k))) && return NaN, NaN
    
    if isnan(segment.k)
        x = segment.b
        y = ray.k * x + ray.b
    elseif isnan(ray.k)
        x = ray.b
        y = segment.k * x + segment.b
    else
        x = (ray.b - segment.b) / (segment.k - ray.k)
        y = ray.k * x + ray.b   
    end
    
    point = [x, y]
    !satisfies_constraints(ray, point) && return NaN, NaN
    !satisfies_constraints(segment, point) && return NaN, NaN
    
    return point, norm(ray.pos - point)
end

function find_intersection(trajectory::Segment, subspace::Segment)
    (isapprox(trajectory.k, subspace.k) || (isnan(trajectory.k) && isnan(subspace.k))) && return NaN
    
    if isnan(trajectory.k)
        x = trajectory.b
        y = subspace.k * x + subspace.b
    elseif isnan(subspace.k)
        x = subspace.b
        y = trajectory.k * x + trajectory.b
    else
        x = (subspace.b - trajectory.b) / (trajectory.k - subspace.k)
        y = subspace.k * x + subspace.b   
    end
    
    point = [x, y]
    ! satisfies_constraints(trajectory, point) && return NaN
    
    point_a = point - subspace.A
    b_a = subspace.B - subspace.A
    t = NaN
    
    if isapprox(b_a[1], 0.0)
        t = point_a[2] / b_a[2]
        ((! isapprox(point[2], subspace.A[2]) && t < 0.0) || (! isapprox(t, 1.0) && t > 1.0)) && return NaN
    else
        t = point_a[1] / b_a[1]
        ((! isapprox(point[1], subspace.A[1]) && t < 0.0) || (! isapprox(t, 1.0) && t > 1.0)) && return NaN
    end
    
    return t
end

function find_intersection(arc::Arc, ray::Ray)
    intersection_points = NaN
    num_of_roots = 0
    
    if isnan(ray.k)
        m, n = arc.C
        num_of_roots, roots = solve_quadratic_equation(1, -2 * n, (ray.b - m)^2 + n^2 - arc.R^2)
        num_of_roots != 0 && (intersection_points = [[ray.b, y] for y ∈ roots])
    else
        m, n = arc.C
        k, b = ray.k, ray.b
        num_of_roots, roots = transition_qe(k, b, m, n, arc.R)
        num_of_roots != 0 && (intersection_points = [[x, k * x + b] for x ∈ roots])
    end
    
    if num_of_roots == 1
        satisfies_constraints(arc, intersection_points[1]) && satisfies_constraints(ray, intersection_points[1]) && return intersection_points[1], norm(ray.pos - intersection_points[1])
    elseif num_of_roots == 2
        dist_first = norm(ray.pos - intersection_points[1]) 
        dist_second = norm(ray.pos - intersection_points[2])
        dist_first > dist_second && begin intersection_points[1], intersection_points[2], dist_first, dist_second = intersection_points[2], intersection_points[1], dist_second, dist_first end
        satisfies_constraints(arc, intersection_points[1]) && satisfies_constraints(ray, intersection_points[1]) && return intersection_points[1], dist_first
        satisfies_constraints(arc, intersection_points[2]) && satisfies_constraints(ray, intersection_points[2]) && return intersection_points[2], dist_second
    end

    return NaN, NaN
end

"""

    calculate_velocity_vector(obstacle, ray, point)

A bunch of functions that calculate the velocity vector for different obstacles
"""
calculate_velocity_vector(segment::Segment, ray::Ray, point::Vector{Float64}) = dot(ray.velocity, segment.s) * segment.s - dot(ray.velocity, segment.n) * segment.n

function calculate_velocity_vector(arc::Arc, ray::Ray, point::Vector{Float64})
    n = (arc.C - point) / norm(arc.C - point) 
    s = isapprox(n[2], 0.0) ? [n[2], -n[1]] : [-n[2], n[1]]
    
    return dot(ray.velocity, s) * s - dot(ray.velocity, n) * n
end

"""

    simulate!(table, path)

Run the simulation on the `table`.
First two rows of `path` represent position, last two velocity.
The first column holds the initial condition, the rest is overwritten. 
"""
function simulate!(table::Table, path::Matrix{Float64})
    num_of_coordinates, num_of_points = size(path)
    (num_of_coordinates != 4 || num_of_points < 2 || isapprox(path[3:4, 1], [0,0])) && throw(ArgumentError("Invalid path argument. Matrix path must have at least two columns and exactly four rows."))

    for point ∈ 2:num_of_points
        ray = Ray(path[1:2, point - 1], path[3:4, point - 1])
        min_dist = Inf
        next_point = NaN
        obstacle_id = 0
    
        for (index, obstacle) ∈ enumerate(table.obstacles)
            intersection_point, dist = find_intersection(obstacle, ray)
            dist < min_dist && begin min_dist, next_point, obstacle_id = dist, intersection_point, index end
        end
        
        min_dist == Inf && throw(ErrorException("RuntimeError: No intersection point")) # There is no RuntimeError in Julia
   
        path[1:2, point], path[3:4, point] = next_point, calculate_velocity_vector(table.obstacles[obstacle_id], ray, next_point)
    end
end

"""

    Base.show(io, x)

Prints a human-readable representation of an obstacle or a table to the specified IO stream.
"""
Base.show(io::IO, x::Segment) = print(io, "Segment: A - ($(x.A[1]), $(x.A[2])), B - ($(x.B[1]), $(x.B[2]))")

Base.show(io::IO, x::Arc) = print(io, "Arc: Center - ($(x.C[1]), $(x.C[2])), radius - $(x.R), α - $(x.α), β - $(x.β)")

Base.show(io::IO, x::Table) = print(io, "Billiard table with $(length(x.obstacles)) obstacles")

"""

    draw_obstacle(obstacle)

Draw obstacles.
"""
function draw_obstacle(segment::Segment)
    move(Luxor.Point(segment.A[1], segment.A[2]))
    line(Luxor.Point(segment.B[1], segment.B[2]))
end

function draw_obstacle(arc_obstacle::Arc)
    α = arc_obstacle.α < 0 ? arc_obstacle.α + 2 * π : arc_obstacle.α
    β = arc_obstacle.β < 0 ? arc_obstacle.β + 2 * π : arc_obstacle.β
    Luxor.arc(Luxor.Point(arc_obstacle.C[1], arc_obstacle.C[2]), arc_obstacle.R, α, β)
end

"""

    draw_path(table, path)

Draw the `table` and a `path`.
"""
function draw_path(table::Table, path::Matrix{Float64})
    @svg begin    
        sethue("red")
        setline(1)
        setlinecap("round")
        setlinejoin("bevel")
        for point ∈ axes(path, 2)
            line(Luxor.Point(path[1, point], path[2, point]))
        end
        strokepath()
        
        sethue("black")
        setline(3)
        for obstacle ∈ table.obstacles
            draw_obstacle(obstacle)
            newsubpath()
        end
        strokepath()
    end
end


"""

    draw_poincare(paths, P1, P2)

Draw the Poincaré section from array of trajectories (`paths`) and defined
by segment with endpoints `P1` and `P2`.

"""
function draw_poincare(paths, P1::Vector{Float64}, P2::Vector{Float64})
    # CairoMakie.activate!()
    fig = Figure(size = (200, 1000)) # it doesn't work, I don't know how to adjust the sizes of figure in another way
    ax = Axis(fig[1, 1])
    hidespines!(ax)
    hidedecorations!(ax)
    ax.title = "Poincarého řezy"
    
    
    subspace = Segment(P1, P2)
    for path ∈ paths
        xs, ys = Float64[], Float64[]
        for point ∈ axes(path, 2)[2:end]
            t = find_intersection(Segment(path[1:2, point - 1], path[1:2, point]), subspace)
            if ! isnan(t)
                velocity = path[3:4, point - 1]
                ϕ = acos(dot(velocity, [0.0, 1.0]) / norm(velocity))
                push!(xs, t)
                push!(ys, ϕ)
            end
        end
        
        scatter!(xs, ys, markersize=3)
    end
    current_figure()
end

end # module
