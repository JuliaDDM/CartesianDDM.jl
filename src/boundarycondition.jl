abstract type BoundaryCondition end

struct Dirichlet <: BoundaryCondition end
struct Neumann <: BoundaryCondition end
struct Periodic <: BoundaryCondition end

