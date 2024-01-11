# Dynamical-Billiards
Theory of chaos studies the behavior of so-called `dynamical billiards`. These are simple models describing the motion of a small ball in a plane confined by obstacles, where perfectly elastic reflections occur upon collision.  In free space, the ball moves uniformly in a straight line, and its direction changes only upon reflection from an obstacle according to the well-known `law of reflection` (the angle of incidence is equal to the angle of reflection). The goal of this project is to create a simulation of these billiards using the Julia programming language.

## Table of contents
- [Theoretical background](#Theoretical-background)
- [Poincare sections](#Poincare-sections)
- [Implementation](#Implementation)

## Theoretical background
In our simulator, we will focus only on two types of obstacles that can be used to construct the edges of the billiard:

* line segments connecting two specified points,
* segments of a circle.

Our ball will be described by its position (a point in the plane) and a velocity vector.
It seems sufficient to record only the points where reflection occurs.

Reflection in both cases occurs similarly. There is only a change in the velocity magnitude $v \in \mathbb{R}^2$ according to the following rule: if $s \in \mathbb{R}^2$ is the vector indicating the tangent direction to the obstacle, and $n \in \mathbb{R}^2$ is the corresponding normal direction, both normalized to unity (length 1), then the new direction of the ball after reflection is

$$
\langle s, v \rangle s - \langle n, v \rangle n,
$$

where $\langle,\rangle$ denotes the standard scalar product. The vectors $s$ and $n$ mentioned above in the case of individual obstacles can be determined easily:

* for a line segment: the direction vector $s$ is the directional vector of the line on which the segment lies, and the normal vector $n$ is any vector perpendicular to $s$ with unit length,
* for a circle segment: the normal vector $n$ is the vector connecting the point of impact and the center of the circle, and the direction vector $s$ is any vector perpendicular to $n$ with unit length.

## Poincare sections
One of the visual tools for observing chaos is the so-called Poincare sections.

The state of the ball in our billiard is determined by its position (a point in the plane) and the direction of velocity (magnitude does not matter). Therefore, any state of the ball is an element of the set $\mathbb{R}^2 \times \langle -\pi, \pi \rangle$ (called the phase space). While elements of this set and trajectories in the plane can be visualized, the result is often quite uncluttered.

The idea behind Poincare is to `"slice"` this set. In our case, we specifically fix a segment (along which we cut the billiard) and only record the intersections of trajectories with this segment and the direction of velocity. If the endpoints of the section are $A$ and $B$, and the trajectory passes through the point $X$ lying on the segment connecting these two points, i.e., $X = A + t (B - A)$, where $t \in \langle 0, 1 \rangle$, then we display a point with coordinates $t, \varphi$, where $\varphi$ is the angle that the velocity vector makes with the positive direction of the $x$-axis, i.e., a value between $-\pi$ and $\pi$.

## Implementation
In the file `src/billiard.jl`, the `Billiard` module is implemented. It's main purpose is to simulate the movement of a ball in a billiard.

The abstract type `Obstacle` models any obstacle and has two subtypes, namely, `Segment` and `Arc`.

* The `Segment` type represents a line segment and requires only two end points of the segment to create an object. The points are considered as `Vector{Float64}` with two components. It is not possible to create a segment between points that are too close or even equal (check the proximity of points using the `isapprox` method).

* The `Arc` type represents a part of a circular arc with a center at a point (`Vector{Float64}`) and a radius (`Float64`), along with two angles (both `Float64`). The first angle must be smaller than the second, indicating the start and end of the arc in a counterclockwise direction from the positive $x$-axis. For example, if we have a unit circle, with angles $0$ and $\pi$, we are referring to only a part above the $x$-axis. The center of the circle must be a two-component vector (as we are in a plane), the radius must be non-negative, and the first angle must be smaller than the second.

The last custom type used is the `Table` type, which straightforwardly models a billiard assembled from obstacles of the two aforementioned types. These obstacles are easily stored/fetched in an array of type `Vector{Obstacle}`.

The `simulate!` function calculates successively the reflection points on the edges of the billiard, which we pass as the first argument (`table`). The second argument (`path`) contains a matrix of `Float64` values of size $4 \times n$, where the first column stores the initial conditions: the coordinates of the position and components of the velocity vector (not necessarily unitary). The method then fills the remaining columns of the array with values of the coordinates of impacts on the boundary and new velocity components (after reflection).

The functions `draw_path` and `draw_poincare` serve for the visualization of dynamic billiards.

The entire implementation, for simplicity and speed, is done in machine numbers (i.e., vectors and various values are considered only of type `Float64`).

This task is created as a project, meaning it can be activated and dependencies used for visualization can be installed (using `Project.toml` and `Manifest.toml`).
