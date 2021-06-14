## code to prepare `SimData` dataset goes here

library(terra)
library(ReacTran)

## =============================================================================
## A 2-D model with diffusion in x- and y direction and first-order
## consumption - more efficient implementation, specifying ALL 2-D grids
## =============================================================================

N     <- 202          # number of grid cells
N2  <- ceiling(N/2)
Dy    <- Dx <- 0.01   # diffusion coeff, X- and Y-direction
ax <- 0.0
ay <- 0.0
r     <- 0.5         # consumption rate
ini   <- 1           # initial value at x=0
k1 = 100.0

x.grid    <- setup.grid.1D(x.up = -5, x.down = 5, N = N)
y.grid    <- setup.grid.1D(x.up = -5, x.down = 5, N = N)
grid2D    <- setup.grid.2D(x.grid, y.grid)

D.grid    <- setup.prop.2D(value = Dx, y.value = Dy, grid = grid2D)
v.grid    <- setup.prop.2D(value = ax, y.value = ay, grid = grid2D)
str(v.grid)
# setting up advection in upper left quadrant
v.grid$x.int[1:102, 1:101] = 0.2
v.grid$y.int[1:101, 1:102] = 0
# setting up advection in upper right quadrant
v.grid$x.int[1:102, 102:202] = 0
v.grid$y.int[1:101, 102:203] = -0.2
# setting up advection in the lower right quadrant
v.grid$x.int[102:203, 102:202] = -0.2
v.grid$y.int[102:202, 102:203] = 0
# setting up advection in the lower left quadrant
v.grid$x.int[102:203, 1:101] = 0
v.grid$y.int[102:202, 1:102] = 0.2

A.grid    <- setup.prop.2D(value = 1, grid = grid2D)
AFDW.grid <- setup.prop.2D(value = 1, grid = grid2D)
VF.grid   <- setup.prop.2D(value = 1, grid = grid2D)

# The model equations - using the grids
DiffAdv2Db <- function (t, y, parms)  {

  CONC  <- matrix(nrow = N, ncol = N, data = y)

  dCONC <- tran.2D(CONC, grid = grid2D, D.grid = D.grid,
                   A.grid = A.grid, VF.grid = VF.grid, AFDW.grid = AFDW.grid,
                   v.grid = v.grid)$dC + r * CONC

  return (list(dCONC))
}

# initial condition: 0 everywhere, except in central point
y <- matrix(nrow = N, ncol = N, data = 0)
# I initialize with a single point in the centre of each
# quadrant.
y[51, 51] <- 1  # initial concentration in the central point...
y[51, 151] <- 1
y[151, 51] <- 1
y[151, 151] <- 1

# solve for 8 time units
times <- 0:8
outb <- ode.2D (y = y, func = DiffAdv2Db, t = times, parms = NULL,
                dim = c(N, N), lrw = 5640000)

# constructing a data-frame
SimData = expand.grid(grid2D$x.mid, grid2D$y.mid)
colnames(SimData) = c("xcoord", "ycoord")
SimData$t1 = as.numeric(subset(outb, subset = (time == 1)))
SimData$t2 = as.numeric(subset(outb, subset = (time == 2)))
SimData$t3 = as.numeric(subset(outb, subset = (time == 3)))
SimData$t4 = as.numeric(subset(outb, subset = (time == 4)))
SimData$t5 = as.numeric(subset(outb, subset = (time == 5)))
SimData$t6 = as.numeric(subset(outb, subset = (time == 6)))

usethis::use_data(SimData, overwrite = TRUE)
