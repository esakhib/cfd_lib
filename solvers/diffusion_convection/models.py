import logging

import numpy as np

from solvers.diffusion_convection.input_data import GridTimeDataDC, InputDataDC
from utils.boundary_type import BoundaryConditionsModel, BoundaryConditionsType
from utils.common import timer
from utils.equation_type import SolverInputData, SolverOutputData
from utils.tdma import run_tdma


class FiniteVolumeModel:
    """Class that contains finite volume method parameters.
    """

    def __init__(self, solver_input_data: SolverInputData):
        # solver input data
        self._solver_input_data: SolverInputData = solver_input_data

        # grid & time scale input data
        self._grid_time_data: GridTimeDataDC = self._solver_input_data.grid_time_data

        self._length: float = self._grid_time_data.x_length
        self._height: float = self._grid_time_data.y_height
        self._total_time: float = self._grid_time_data.total_time

        self._nx: int = self._grid_time_data.nx
        self._ny: int = self._grid_time_data.ny
        self._nt: int = self._grid_time_data.nt

        self._dx: float = self._length / (self._nx - 1)
        self._dy: float = self._height / (self._ny - 1)

        # unsteady problem
        if self._nt - 1 != 0:
            self._dt: float = self._total_time / (self._nt - 1)
        # steady-state problem
        else:
            self._dt: float = 0.0

        self._dx_e: np.ndarray = np.full(
            shape=(self._nx, self._ny),
            fill_value=self._dx,
            dtype=np.float64
        )

        self._dx_w: np.ndarray = np.full(
            shape=(self._nx, self._ny),
            fill_value=self._dx,
            dtype=np.float64
        )

        self._dy_n: np.ndarray = np.full(
            shape=(self._nx, self._ny),
            fill_value=self._dy,
            dtype=np.float64
        )

        self._dy_s: np.ndarray = np.full(
            shape=(self._nx, self._ny),
            fill_value=self._dy,
            dtype=np.float64
        )

    @property
    def length(self) -> float:
        return self._length

    @property
    def height(self) -> float:
        return self._height

    @property
    def total_time(self) -> float:
        return self._total_time

    @property
    def nx(self) -> int:
        return self._nx

    @property
    def ny(self) -> int:
        return self._ny

    @property
    def nt(self) -> int:
        return self._nt

    @property
    def dx(self) -> float:
        return self._dx

    @property
    def dy(self) -> float:
        return self._dy

    @property
    def dt(self) -> float:
        return self._dt

    @property
    def dx_e(self) -> np.ndarray:
        return self._dx_e

    @property
    def dx_w(self) -> np.ndarray:
        return self._dx_w

    @property
    def dy_n(self) -> np.ndarray:
        return self._dy_n

    @property
    def dy_s(self) -> np.ndarray:
        return self._dy_s


class PhysicsModel:
    """Class that contains physics equation parameters.
    """

    def __init__(self, solver_input_data: SolverInputData, fvm_model: FiniteVolumeModel):
        # input data
        self._solver_input_data: SolverInputData = solver_input_data

        # finite volume method model
        self._fvm_model: FiniteVolumeModel = fvm_model

        # equation input data
        self._equation_input_data: InputDataDC = self._solver_input_data.equation_input_data

        # boundary conditions model
        self._boundary_conditions_model: BoundaryConditionsModel = self._equation_input_data.boundary_conditions_model

        # initial value
        self._c_init: float = self._equation_input_data.c_init

        # boundary values
        self._c_left: float = self._equation_input_data.c_left
        self._c_right: float = self._equation_input_data.c_right
        self._c_upper: float = self._equation_input_data.c_upper
        self._c_lower: float = self._equation_input_data.c_lower

        # source value
        self._q_source: float = self._equation_input_data.q_source

        # diffusion coefficient
        self._d_coef: float = self._equation_input_data.d

        self._d_e: np.ndarray = np.full(
            shape=(self._fvm_model.nx, self._fvm_model.ny),
            fill_value=self._d_coef,
            dtype=np.float64
        )

        self._d_w: np.ndarray = np.full(
            shape=(self._fvm_model.nx, self._fvm_model.ny),
            fill_value=self._d_coef,
            dtype=np.float64
        )

        self._d_n: np.ndarray = np.full(
            shape=(self._fvm_model.nx, self._fvm_model.ny),
            fill_value=self._d_coef,
            dtype=np.float64
        )

        self._d_s: np.ndarray = np.full(
            shape=(self._fvm_model.nx, self._fvm_model.ny),
            fill_value=self._d_coef,
            dtype=np.float64
        )

        self._u_sed_e, self._u_sed_w, self._u_sed_n, self._u_sed_s = [None] * 5

        if self._equation_input_data.use_velocity:
            self._u_sed_e: np.ndarray = np.zeros(
                shape=(self._fvm_model.nx, self._fvm_model.ny),
                dtype=np.float64
            )

            self._u_sed_w: np.ndarray = np.zeros(
                shape=(self._fvm_model.nx, self._fvm_model.ny),
                dtype=np.float64
            )

            self._u_sed_s: np.ndarray = np.zeros(
                shape=(self._fvm_model.nx, self._fvm_model.ny),
                dtype=np.float64
            )

            self._u_sed_n: np.ndarray = np.zeros(
                shape=(self._fvm_model.nx, self._fvm_model.ny),
                dtype=np.float64
            )

    def calc_u_sed(self, c: np.ndarray) -> np.ndarray:
        """Calculate Used by C values.
        """

        u_sed: np.ndarray = np.zeros(
            shape=(self._fvm_model.nx, self._fvm_model.ny),
            dtype=np.float64
        )

        const = self._equation_input_data.const_u_sed
        r0 = self._equation_input_data.r0
        g = self._equation_input_data.g
        rho1 = self._equation_input_data.rho1
        rho2 = self._equation_input_data.rho2
        mu2 = self._equation_input_data.mu2
        f_c = self._equation_input_data.f_c

        for i in np.arange(1, self._fvm_model.nx - 1):
            u_sed[i] = const * r0 ** 2.0 * g * (rho1 - rho2) * f_c((c[i] + c[i + 1]) / 2.0) / mu2

        return u_sed

    def update_u_sed(self, u_sed: np.ndarray) -> None:
        """Update velocity.
        """

        self._u_sed_w = u_sed
        self._u_sed_e = u_sed
        self._u_sed_n = u_sed
        self._u_sed_s = u_sed

    @property
    def c_init(self) -> float:
        return self._c_init

    @property
    def c_left(self) -> float:
        return self._c_left

    @property
    def c_right(self) -> float:
        return self._c_right

    @property
    def c_upper(self) -> float:
        return self._c_upper

    @property
    def c_lower(self) -> float:
        return self._c_lower

    @property
    def q_source(self) -> float:
        return self._q_source

    @property
    def boundary_conditions_model(self) -> BoundaryConditionsModel:
        return self._boundary_conditions_model

    @property
    def d_e(self) -> np.ndarray:
        return self._d_e

    @property
    def d_w(self) -> np.ndarray:
        return self._d_w

    @property
    def d_s(self) -> np.ndarray:
        return self._d_s

    @property
    def d_n(self) -> np.ndarray:
        return self._d_n

    @property
    def u_sed(self) -> np.ndarray:
        return self._u_sed

    @property
    def u_sed_e(self) -> np.ndarray:
        return self._u_sed_e

    @property
    def u_sed_w(self) -> np.ndarray:
        return self._u_sed_w

    @property
    def u_sed_n(self) -> np.ndarray:
        return self._u_sed_n

    @property
    def u_sed_s(self) -> np.ndarray:
        return self._u_sed_s


class DiscreteAnalogueModel:
    """Class that contains discrete analouge scheme.
    """

    def __init__(self, fvm_model: FiniteVolumeModel, physics_model: PhysicsModel):
        self._fvm_model: FiniteVolumeModel = fvm_model
        self._physics_model: PhysicsModel = physics_model

        # left-hand side
        self._a_p: np.ndarray = np.zeros(
            shape=(self._fvm_model.nx, self._fvm_model.ny),
            dtype=np.float64
        )

        self._a_e: np.ndarray = np.zeros(
            shape=(self._fvm_model.nx, self._fvm_model.ny),
            dtype=np.float64
        )

        self._a_w: np.ndarray = np.zeros(
            shape=(self._fvm_model.nx, self._fvm_model.ny),
            dtype=np.float64
        )

        self._a_n: np.ndarray = np.zeros(
            shape=(self._fvm_model.nx, self._fvm_model.ny),
            dtype=np.float64
        )

        self._a_s: np.ndarray = np.zeros(
            shape=(self._fvm_model.nx, self._fvm_model.ny),
            dtype=np.float64
        )

        # right-hand side
        self._b: np.ndarray = np.zeros(
            shape=(self._fvm_model.nx, self._fvm_model.ny),
            dtype=np.float64
        )

    def apply_boundary_conditions(self):
        c_left: float = self._physics_model.c_left
        c_right: float = self._physics_model.c_right
        q_source: float = self._physics_model.q_source
        d_e: np.ndarray = self._physics_model.d_e
        dx_e: np.ndarray = self._fvm_model.dx_e
        d_w: np.ndarray = self._physics_model.d_w
        dx_w: np.ndarray = self._fvm_model.dx_w
        u_sed_e: np.ndarray = self._physics_model.u_sed_e
        u_sed_w: np.ndarray = self._physics_model.u_sed_w

        if self._physics_model.boundary_conditions_model.left_side == BoundaryConditionsType.Dirichlet:
            self._a_e[0] = -1.0
            self._a_w[0] = 0.0
            self._a_p[0] = 1.0
            self._b[0] = 2.0 * c_left

        if self._physics_model.boundary_conditions_model.left_side == BoundaryConditionsType.Neumann:
            self._a_e[0] = 1.0
            self._a_w[0] = 0.0
            self._a_p[0] = 1.0
            self._b[0] = -q_source / d_e[0] * dx_e[0]

        if self._physics_model.boundary_conditions_model.left_side == BoundaryConditionsType.Robin:
            self._a_e[0] = 1.0 - u_sed_e[0] * dx_e[0] / (2.0 * d_e[0])
            self._a_w[0] = 0.0
            self._a_p[0] = 1.0 + u_sed_e[0] * dx_e[0] / (2.0 * d_e[0])
            self._b[0] = 0.0

        if self._physics_model.boundary_conditions_model.right_side == BoundaryConditionsType.Dirichlet:
            self._a_e[-1] = 0.0
            self._a_w[-1] = -1.0
            self._a_p[-1] = 1.0
            self._b[-1] = 2.0 * c_right

        if self._physics_model.boundary_conditions_model.right_side == BoundaryConditionsType.Neumann:
            self._a_e[-1] = 0.0
            self._a_w[-1] = 1.0
            self._a_p[-1] = 1.0
            self._b[-1] = q_source / d_w[-1] * dx_w[-1]

        if self._physics_model.boundary_conditions_model.right_side == BoundaryConditionsType.Robin:
            self._a_e[-1] = 0.0
            self._a_w[-1] = 1.0 + u_sed_w[-1] * dx_w[-1] / (2.0 * d_w[-1])
            self._a_p[-1] = 1.0 - u_sed_w[-1] * dx_w[-1] / (2.0 * d_w[-1])
            self._b[-1] = 0.0

    def make_discrete_analogue(self, solution_data: np.ndarray):
        for i in np.arange(1, self._fvm_model.nx - 1):
            self._a_e[i] = (self._physics_model.d_e[i] / self._fvm_model.dx_e[i] +
                            max(-self._physics_model.u_sed_e[i], 0.0))

            self._a_w[i] = (self._physics_model.d_w[i] / self._fvm_model.dx_w[i] +
                            max(self._physics_model.u_sed_w[i], 0.0))

            self._a_p[i] = (self._fvm_model.dx / self._fvm_model.dt +
                            self._physics_model.d_e[i] / self._fvm_model.dx_e[i] +
                            self._physics_model.d_w[i] / self._fvm_model.dx_w[i] +
                            (self._physics_model.u_sed_e[i] - self._physics_model.u_sed_w[i]))

            self._b[i] = self._fvm_model.dx / self._fvm_model.dt * solution_data[i]

    @property
    def a_p(self) -> np.ndarray:
        return self._a_p

    @a_p.setter
    def a_p(self, _a_p: np.ndarray):
        self._a_p = _a_p

    @property
    def a_e(self) -> np.ndarray:
        return self._a_e

    @a_e.setter
    def a_e(self, _a_e: np.ndarray):
        self._a_e = _a_e

    @property
    def a_w(self) -> np.ndarray:
        return self._a_w

    @a_w.setter
    def a_w(self, _a_w: np.ndarray):
        self._a_w = _a_w

    @property
    def a_n(self) -> np.ndarray:
        return self._a_n

    @a_n.setter
    def a_n(self, _a_n: np.ndarray):
        self._a_n = _a_n

    @property
    def a_s(self) -> np.ndarray:
        return self._a_s

    @a_s.setter
    def a_s(self, _a_s: np.ndarray):
        self._a_s = _a_s

    @property
    def b(self) -> np.ndarray:
        return self._b

    @b.setter
    def b(self, _b: np.ndarray):
        self._b = _b


class DiffusionConvectionSolverModel:
    """Solver for diffusion-convection problem.
    """

    @timer
    def __init__(self, solver_input_data: SolverInputData):
        logging.info('Start initialization grid and solver data ...')

        # input data
        self._solver_input_data: SolverInputData = solver_input_data

        # output data
        self._solver_output_data: SolverOutputData = SolverOutputData()

        # finite volume method model
        self._fvm_model: FiniteVolumeModel = FiniteVolumeModel(
            solver_input_data=self._solver_input_data
        )

        # physics model
        self._physics_model: PhysicsModel = PhysicsModel(
            solver_input_data=self._solver_input_data,
            fvm_model=self._fvm_model
        )

        # discrete analogue model
        self._discrete_analogue_model: DiscreteAnalogueModel = DiscreteAnalogueModel(
            fvm_model=self._fvm_model,
            physics_model=self._physics_model
        )

        # steady-state or unsteady problem
        self._is_unsteady = True if self._fvm_model.dt == 0 else False

        self._current_time_step_solution: np.ndarray = np.zeros(
            shape=(self._fvm_model.nx, self._fvm_model.ny),
            dtype=np.float64
        )

        self._old_time_step_solution: np.ndarray = np.zeros(
            shape=(self._fvm_model.nx, self._fvm_model.ny),
            dtype=np.float64
        )

        # solution data by each time iteration
        self._solutions: dict = {}

        # velocity data by each time iteration
        self._velocity: dict = {}

        # domain output grid
        self._solver_output_data.grid = np.arange(
            start=0.0,
            stop=self._fvm_model.length,
            step=self._fvm_model.dx
        )

        # add Length value at domain grid end
        self._solver_output_data.grid = np.append(
            arr=self._solver_output_data.grid,
            values=self._fvm_model.length
        )

        self._solver_output_data.numerical_solution = self._current_time_step_solution
        self._solver_output_data.total_solutions = self._solutions
        self._solver_output_data.total_velocity = self._velocity

        logging.info('End initialization grid and solver data.')

    def get_next_time_step_solution(self):
        """Solve the equation by TDMA algorithm.
        """

        return run_tdma(
            a=self._discrete_analogue_model.a_p,
            b=self._discrete_analogue_model.a_e,
            c=self._discrete_analogue_model.a_w,
            d=self._discrete_analogue_model.b
        )

    @timer
    def solve_numerical(self):
        """Return numerical solution from FiniteVolumeScheme.
        """

        if self._is_unsteady:
            logging.info('Start numerical solution ...')

            # set initial condition
            self._old_time_step_solution = np.full_like(
                a=self._current_time_step_solution,
                fill_value=self._physics_model.c_init
            )

            # start time, sec
            current_time = 0.0

            # loop through all time layers
            while current_time <= self._fvm_model.total_time:
                # use velocity
                if self._solver_input_data.equation_input_data.use_velocity:
                    # calculate velocity
                    u_sed = self._physics_model.calc_u_sed(
                        c=self._old_time_step_solution
                    )

                    # update velocity
                    self._physics_model.update_u_sed(
                        u_sed=u_sed
                    )

                # initialize discrete analogue scheme using data at the current time layer
                self._discrete_analogue_model.make_discrete_analogue(
                    solution_data=self._old_time_step_solution
                )

                # apply boundary conditions
                self._discrete_analogue_model.apply_boundary_conditions()

                # get data at the next time layers
                self._current_time_step_solution = self.get_next_time_step_solution()

                # update solution
                self._old_time_step_solution = self._current_time_step_solution

                # store data in dict
                self._solutions[current_time] = self._current_time_step_solution
                self._velocity[current_time] = self._physics_model.u_sed

                # go to the next time layer
                current_time += self._fvm_model.dt

        logging.info('End numerical solution.')

    @property
    def current_time_step_solution(self) -> np.ndarray:
        return self.current_time_step_solution

    @property
    def solver_output_data(self) -> SolverOutputData:
        return self._solver_output_data
