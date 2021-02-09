import sympy
from julia import Pkg
Pkg.activate(".\\JuliaFunctions")
from julia import JuliaFunctions

print("done with julia")


class VaporSorptionSystem:
    def __init__(self, pressure_dependent_error=True):

        self.pressure_dependent_error = pressure_dependent_error

        # define known variables
        self.vc = sympy.Symbol("V_c")       # charge chamber volume
        self.vs = sympy.Symbol("V_s")       # sampling chamber volume
        self.R = sympy.Symbol("R")
        self.T = sympy.Symbol("T_K")
        self.pen_mw = sympy.Symbol("MW_Pen")
        self.pol_dry_mass = sympy.Symbol("(\"pol. dry mass\")")
        self.pol_dens = sympy.Symbol("rho_pol")

        self.n_p_variable = sympy.Symbol("n_p")        # moles absorbed by the polymer
        self.n_p_minus1 = sympy.Symbol("n_(p-1)")  # moles absorbed by the polymer in the previous step

        # pressure measurements
        self.pi = sympy.Symbol("P_i")  # initial pressure of system
        self.pf = sympy.Symbol("P_f")  # final pressure of system
        self.pf_minus1 = sympy.Symbol("P_(f-1)")  # final system pressure in previous step
        self.pcf_minus1 = sympy.Symbol("P_(cf-1)")  # charge chamber pressure after valve is closed in previous step

        # measured / found variances
        if self.pressure_dependent_error:
            self.sigma_p = sympy.Symbol("\\sigma_{%P}")
            self.sigma_p_pi = self.sigma_p * self.pi
            self.sigma_p_pf = self.sigma_p * self.pf
            self.sigma_p_pf_minus1 = self.sigma_p * self.pf_minus1
            self.sigma_p_pcf_minus1 = self.sigma_p * self.pcf_minus1
        else:  # else: if the pressure variance is a static number, then all measurements from it have the same error
            self.sigma_p = sympy.Symbol("\\sigma_{P}")
            self.sigma_p_pi = self.sigma_p
            self.sigma_p_pf = self.sigma_p
            self.sigma_p_pf_minus1 = self.sigma_p
            self.sigma_p_pcf_minus1 = self.sigma_p

        self.sigma_v_c = sympy.Symbol("sigma_v_c")
        self.sigma_v_s = sympy.Symbol("sigma_v_s")
        self.sigma_T = sympy.Symbol("sigma_T")
        self.sigma_np = sympy.Symbol("sigma_np")
        self.sigma_np_minus1 = sympy.Symbol("sigma_(np-1)")
        self.sigma_pol_dry_mass = sympy.Symbol("sigma_(\"pol. dry mass\")")
        self.sigma_pol_density = sympy.Symbol("\\sigma_{\\rho_{pol}}")

        # variable-variance pairs
        self.concentration_variable_dict = {
            self.n_p_variable: self.sigma_np,
            self.pol_dry_mass: self.sigma_pol_dry_mass,
            self.pol_dens: self.sigma_pol_density,
        }
        self.np_variable_dict = {
            self.pi: self.sigma_p_pi,
            self.pf: self.sigma_p_pf,
            self.pf_minus1: self.sigma_p_pf_minus1,
            self.pcf_minus1: self.sigma_p_pcf_minus1,
            self.T: self.sigma_T,
            self.vc: self.sigma_v_c,
            self.vs: self.sigma_v_s,
            self.n_p_minus1: self.sigma_np_minus1
        }

        # found expressions
        self.concentration = self._concentration()

    def _init_mol_in_charge_chamber(self):
        return self.pi * (self.vc*(0.01**3)) / (self.R * self.T)

    def _final_mol_in_system(self):
        return self.pf * ((self.vc + self.vs) * (0.01 ** 3)) / (self.R * self.T)

    def _final_mol_in_system_at_prev_step(self):
        return self.pf_minus1 * ((self.vc + self.vs) * (0.01 ** 3)) / (self.R * self.T)

    def _final_mol_in_charge_chamber_at_prev_step(self):
        return self.pcf_minus1 * (self.vc * (0.01 ** 3)) / (self.R * self.T)

    def n_p(self):
        n_c0 = self._init_mol_in_charge_chamber()
        n_f = self._final_mol_in_system()
        n_f_minus1 = self._final_mol_in_system_at_prev_step()
        n_c_minus1 = self._final_mol_in_charge_chamber_at_prev_step()
        n_p_minus1 = self.n_p_minus1
        return n_c0 - n_f + n_f_minus1 - n_c_minus1 + n_p_minus1

    def _concentration(self):
        return self.n_p_variable * self.pen_mw / self.pol_dry_mass * 22414 * self.pol_dens / self.pen_mw

    def concentration_variance(self):
        concentration_variance = sympy.simplify(ErrorPropagation.determine_variance(
            self.concentration,
            self.concentration_variable_dict))
        print(sympy.latex(concentration_variance))

    def n_p_variance(self):
        if self.pressure_dependent_error:
            print("Pressure dependent n_p error: ", end='')
        else:
            print("Pressure independent n_p error: ", end='')
        print()

        n_p_variance = ErrorPropagation.determine_variance(
            self.n_p(),
            self.np_variable_dict)
        print("Simplified:", sympy.latex(sympy.simplify(n_p_variance)))
        print("Raw:", sympy.latex(n_p_variance))
        print()


class GasSorptionSystem:
    class Constants:
        B = sympy.Symbol("B")
        A = sympy.Symbol("A")
        R = sympy.Symbol("R")
        T = sympy.Symbol("T")

    def __init__(self):
        pass
        # self.vars = Vars()

    class Vars:
        # just a convenient wrapper for all the variable definitions
        def __init__(self):
            pass

    class PengRobinsonEquation:
        def __init__(self, var_subscript, constants):
            self.vmol = sympy.Symbol("vmol_" + var_subscript)
            self.pressure = sympy.Symbol("p_" + var_subscript)
            self.constants = constants

        def original_objective_function(self):
            """
                Original objective function as seen in the PR expression within the gas sorption template
                We need to solve this analytically for molar volume of the system
                In sympy, if you try to solve an expression, it implicitly sets the expression to zero to solve

            """
            # from VBA code:
            # expression = vmol * (
            #   pressure * vmol ^ 2
            #   + (pressure * B - R * temperature) * vmol
            #   + (A - 2 * R * temperature * B - 3 * pressure * (B ^ 2))
            # )
            # + pressure * B ^ 3 + R * temperature * B ^ 2 - A * B

            objective = self.vmol * (
                    self.pressure * self.vmol**2
                    + (self.pressure * self.constants.B - self.constants.R * self.constants.T) * self.vmol
                    + (self.constants.A - 2 * self.constants.R * self.constants.T * self.constants.B - 3 *
                       self.pressure * (self.constants.B**2))
                    )\
                + self.pressure * self.constants.B ** 3 + self.constants.R * self.constants.T * self.constants.B ** 2 \
                - self.constants.A * self.constants.B
            return objective

        def get_expression_for_moles(self):
            print(sympy.latex(self.original_objective_function()))
            result = sympy.solve(self.original_objective_function(), self.vmol)
            print(len(result))
            print(sympy.latex(sympy.simplify(result[0])))

    def test(self):
        self.PengRobinsonEquation("test", self.Constants).get_expression_for_moles()


class ErrorPropagation:

    @staticmethod
    def determine_variance(equation, variable_variance_dict):
        """
        :param equation: Sympy expression (not an equals() expression, rather, the "something" in f(x)=something
        :param variable_variance_dict: dict of sympy symbols or expressions
                formatted as variables and their corresponding variances, i.e., {variable: variance, ...}
                All variables referenced MUST be present as the same object in memory, or they will not be recognized.
        :type variable_variance_dict: dict
        :return: Sympy expression of the propagated error
        """
        variable_derivative_dict = {}
        for variable in variable_variance_dict.keys():
            variable_derivative_dict[variable] = sympy.diff(equation, variable)

        dummy = sympy.Symbol("dummy")  # I honestly, cannot remember the easier way to create an empty expression
        result = dummy

        for variable in variable_derivative_dict.keys():
            result = result + (variable_variance_dict[variable] ** 2) * (variable_derivative_dict[variable]**2)

        result = result - dummy  # an equally stupid way for me to do this
        result = sympy.sqrt(result)
        return result


GasSorptionSystem().test()  # fundamentally should not work
# VaporSorptionSystem().concentration_variance()
# VaporSorptionSystem().n_p_variance()
# VaporSorptionSystem(pressure_dependent_error=False).n_p_variance()