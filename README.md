# IFGF-RP
## Construction and Initialization Pipeline

The `Solver` class is configured in **three stages**:

1.  **Construct the solver object** using one of two constructors. This save the geometry type.
2.  **Optionally set parameters** 
3.  **Call `init_solver()`**, which finalizes MPI, geometry,
    discretization, and internal data structures.

------------------------------------------------------------------------

## 1. Constructors

The solver provides two construction modes, depending on whether the
geometry is spherical or loaded from a directory.

------------------------------------------------------------------------

## **Constructor 1 --- Sphere Geometry**

``` cpp
template<int PS, int PT>
Solver<PS,PT>::Solver(
    double sphere_radius,
    double sphere_centerX, double sphere_centerY, double sphere_centerZ);
```

### Purpose

Constructs a built-in spherical geometry.

### Behavior

-   Sets:

    ``` cpp
    SPHERE_RADIUS = sphere_radius;
    SPHERE_CENTER = {sphere_centerX, sphere_centerY, sphere_centerZ};
    ```
------------------------------------------------------------------------

## **Constructor 2 --- Geometry Loaded From Directory**

``` cpp
template<int PS, int PT>
Solver<PS,PT>::Solver(const std::string directory,
                      const std::string file_prefix);
```

### Purpose

Constructs the solver by loading geometry from files located in a
directory. Each file corresponds to one geometry patch. These must be stored as director + file_prefix + i where i is the number ranging over the patches. 

## 2. Setting Parameters (Optional)

After constructing the solver, you may configure various parameters (see below):

These settings may be applied in any order **before calling
`init_solver()`**.

------------------------------------------------------------------------

## 3. Initialization Routine

``` cpp
template<int PS, int PT>
void Solver<PS,PT>::init_solver(
    const bool timing,
    const std::complex<double> k,
    MPI_Comm mpi_comm);
```

## 4. Full Usage Example

``` cpp
// 1. Construct solver (sphere example)
Solver solver(1.0, 0.0, 0.0, 0.0);

// 2. Configure solver
solver.set_eq_form(3);                     // Combined layer
solver.set_int_ext(1);                     // Exterior problem
solver.set_n_pts_per_patch(8, 8);          // Increase resolution
solver.set_num_split_per_patch(20, 20);    // Refine patch partitioning
solver.set_coup_param({0.1, 0.0});         // Coupling parameter

// 3. Initialize
solver.init_solver(
    /*timing=*/true,
    /*k=*/std::complex<double>(M_PI,0),
    /*mpi_comm=*/MPI_COMM_WORLD
);
```


## Class Parameter Guide

Below we describe all user configurable values 

## MPI and High-Level Info

-   `MPI_Comm get_mpi_comm() const`: Returns the MPI communicator.
-   `long long get_num_unknowns() const`: Total number of unknowns.
-   `int get_patch_split_x() const`, `get_patch_split_y() const`: Patch
    partition counts.
-   `int get_nlevels_IFGF() const`: Number of IFGF levels.
-   `std::complex<double> get_coup_param() const`: Coupling parameter.

## Interpolation Order
- `Solver<PS,PT>` This is set via the template parameters, where PS is the number of distance interpolation 
points, and PT is the number of angular interpolation points

## Integral Equation Formulation

`set_eq_form(int form)`
Values:
1 = Single Layer
2 = Double Layer
3 = Combined Layer (default)
4 = Other

## Interior/Exterior Problem

`set_int_ext(int int_ext)`
- -1 = Interior
- 1 = Exterior (default)

## Discretization Parameters
-  
-   `set_n_pts_per_patch(p1, p2)` sets points per patch. Default
    `{6,6}`
-   `set_num_wl_per_patch(double)` sets wavelength-based subdivision. This is used by default with 1 wave length per patch direction on each side for non-sphereical geometries.
-   `set_num_split_per_patch(s1, s2)` sets patch subdivisions, and doesn't use wavelength-based subdivision. Default 
    `{16,16}`

-   `get_disc_points_x/y/z()` return global discretization nodes.

## Singular / Near-Singular Integration

-   `N_PTS_SING_INT` Number of points used in the singular/near singular integration. Default `{40,40}`
-   `DELTA_METHOD` Use either proximity box size `1` or percentage box size `2` to determine near singular integrals. Default `2` 
-   `PROXIMITY_BOX_SIZE = 0.1/16` Use a specified box size to determine near singular points
-   `PERCENT_BOX_SIZE = 0.15` Use a percentage of the box size to determine the near singular points.

## GMRES Solver Options

-   `MAX_ITER = 100`
-   `TOL_GMRES = 1e-4`

## Geometry Options

-   `GEOMETRY = 0` → sphere
-   Sphere parameters: `SPHERE_RADIUS`, `SPHERE_CENTER`
-   Edge flags exist but are currently unused.

## Incident Field Parameters

`PLANE_OR_POINT = 0`
- 0 = plane wave
- 1 = point sources

Plane wave: `WAVE_NUMBER`, `PLANE_WAVE_THE`, `PLANE_WAVE_PHI`
- The plane wave travels at 

Point sources: 
- `NUM_POINT_SOURCES` The number of point sources used in the computation.

- `POINT_SOURCE_CENTER` A `std::vector<std::vector<double>>` which holds x,y,z locations of the point sources.

## IFGF Acceleration

-   `USE_ACCELERATOR = true`
-   `N_LEVELS_IFGF = 5` Not used if adaptivity is used.
-   `USE_ADAPTIVITY = false`
-   `MAX_ELEMS_LEAF = 100`

## OpenMP Options

-   `NTHREADS = 16`

## Other Setters
- `set_coup_param(std::complex<double> cp)` Sets the complex coupling parameter. The coupling parameter is chosen automatically based on the wavenumber otherwise.
- `set_num_split_per_patch(double split_1, double split_2)` Sets the number of patch splits. Default is to do this automatically based on the wavenumber for non sphere geometries, or defaults to `{16,16}` for the sphere.
- `set_n_pts_per_patch(double pts1, double pts2)` Sets the number of points to use in each patch. Default `{6,6}`


