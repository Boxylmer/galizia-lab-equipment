from ErrorPropagation import ErrorPropagation
import sympy
import numba as nb
import matplotlib.pyplot as plt

# todo swap out variance for uncertainty


def plotstuff(x1, y1, data2):
    fig, axes = plt.subplots(2, constrained_layout=True)
    # fig.tight_layout()
    fig.suptitle(str(VaporSorptionSystem.DEFAULT_ITERATIONS) + ' Monte-Carlo Iterations', fontsize=12)
    axes[0].scatter(x1, y1, alpha=0.002, s=1)
    axes[0].set_xlabel('Pressure (Pascal)')
    axes[0].set_ylabel('Conc (cc(STP)/cc)')
    axes[0].set_xlim(2600, 2750)
    axes[0].set_ylim(29, 31)

    axes[1].hist(data2, bins=100, orientation='vertical')
    axes[1].set_xlabel('Concentration histogram')

    plt.show()


# first example: deterministic (use ErrorPropagation.analytical)
#                in this case, we may use normal error propagation with symbolic math (SymPy)
class _VaporSorptionEquations:
    def __init__(self, pressure_dependent_error):

        self.pressure_dependent_error = pressure_dependent_error

        # define known variables
        self.vc = sympy.Symbol("V_c")       # charge chamber volume (SPECIFIED AS CM3)
        self.vs = sympy.Symbol("V_s")       # sampling chamber volume
        self.R = sympy.Symbol("R")          # (SPECIFIED AS J/MolK or LBar/MolK)
        self.T = sympy.Symbol("T_K")
        self.pen_mw = sympy.Symbol("MW_Pen")
        self.pol_dry_mass = sympy.Symbol("(\"pol. dry mass\")")
        self.pol_dens = sympy.Symbol("rho_pol")

        self.n_p_variable = sympy.Symbol("n_p")    # moles absorbed by the polymer
        self.n_p_minus1 = sympy.Symbol("n_(p-1)")  # moles absorbed by the polymer in the previous step

        # pressure measurements (ALL ARE PASCAL INPUTS)
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
        concentration_variance = sympy.simplify(ErrorPropagation.create_analytical_error_function(
            self.concentration,
            self.concentration_variable_dict))
        return sympy.latex(concentration_variance)

    def n_p_variance(self):
        n_p_variance = ErrorPropagation.create_analytical_error_function(
            self.n_p(),
            self.np_variable_dict)
        return sympy.latex(sympy.simplify(n_p_variance))

    def get_concentration_function(self):
        return sympy.lambdify([self.n_p_variable, self.pen_mw, self.pol_dry_mass, self.pol_dens],
                              self.concentration)

    def get_n_p_function(self):
        return sympy.lambdify(
            [self.R, self.T, self.pi, self.pf, self.vc, self.vs, self.pf_minus1, self.pcf_minus1, self.n_p_minus1],
            self.n_p())


class VaporSorptionSystem:
    # define some constants
    R = 8.314
    VC = 29.4778313316133
    VS = 7.613879765219

    PRESSURE_MEASUREMENT_VARIANCE = 0.15 / 100  # percent -> multiplier
    TEMPERATURE_MEASUREMENT_VARIANCE = 1  # K
    VC_VARIANCE = 0.098002  # cm3
    VS_VARIANCE = 0.023176  # cm3

    DEFAULT_ITERATIONS = 100000
    USE_NORMAL_DISTRIBUTION = True

    def __init__(self,
                 polymer_density, polymer_mass, temperature, penetrant_mw,
                 polymer_density_error, polymer_mass_error, temperature_error, pen_mw_error):

        self.analytical_np_function = self._get_analytical_np_function()
        self.analytical_concentration_function = self._get_analytical_concentration_function()

        self.polymer_density, self.polymer_density_error = polymer_density, polymer_density_error  # g/cm3
        self.polymer_mass, self.polymer_mass_error = polymer_mass, polymer_mass_error  # g
        self.temperature, self.temperature_error = temperature, temperature_error  # K
        self.pen_mw, self.pen_mw_error = penetrant_mw, pen_mw_error

        self.initial_pressures = []  # pa
        self.final_pressures = []    # pa
        self.final_charge_chamber_pressures = []  # pa

        self.calculated_nps = []
        self.calculated_concentrations = []
        self.calculated_np_errors = []
        self.calculated_concentration_errors = []

    def add_step_data(self, pi, pf, pcf):
        self.initial_pressures.append(pi)
        self.final_pressures.append(pf)
        self.final_charge_chamber_pressures.append(pcf)

    def get_monte_carlo_isotherm(self):
        # todo assert used variables are actually specified
        assert len(self.initial_pressures) == len(self.final_pressures)
        assert len(self.initial_pressures) == len(self.final_charge_chamber_pressures)

        n_steps = len(self.initial_pressures)
        for step_idx in range(n_steps):
            print("Starting step " + str(step_idx))
            pi = self.initial_pressures[step_idx]
            pf = self.final_pressures[step_idx]
            # pcf = self.final_charge_chamber_pressures[step_idx]  # we don't actually need this for the immediate step

            if step_idx == 0:  # if we're evaluating the first step in the isotherm
                pf_minus1 = 0
                pcf_minus1 = 0
                np_minus_1 = 0
                np_minus_1_uncertainty = 0
            else:
                pf_minus1 = self.final_pressures[step_idx - 1]
                pcf_minus1 = self.final_charge_chamber_pressures[step_idx - 1]
                np_minus_1 = self.calculated_nps[step_idx - 1]
                np_minus_1_uncertainty = self.calculated_np_errors[step_idx - 1]

            (np_analytical, np_uncertainty, np_mc_outputs), \
                (conc_analytical, conc_uncertainty, conc_mc_outputs) = self.get_monte_carlo_step(
                pi=pi, pf=pf, pf_minus1=pf_minus1, pcf_minus1=pcf_minus1,
                np_minus_1=np_minus_1, np_minus_1_uncertainty=np_minus_1_uncertainty)
            self.calculated_nps.append(np_analytical)
            self.calculated_np_errors.append(np_uncertainty)
            self.calculated_concentrations.append(conc_analytical)
            self.calculated_concentration_errors.append(conc_uncertainty)
            print("STEP", step_idx, "results", "\n", "Pi:", pi, "\n", "Pf:", pf, "\n",  # "Pcf:", pcf, "\n",
                  "pf_minus1:", pf_minus1, "\n", "pcf_minus1:", pcf_minus1, "\n", "np_minus_1:", np_minus_1, "\n",
                  "np_minus_1_uncertainty:", np_minus_1_uncertainty, "\n",

                  "Calculated Concentration:", self.calculated_concentrations[step_idx], "\n",
                  "Concentration uncertainty:", self.calculated_concentration_errors[step_idx], "\n",
                  "Calculated np:", self.calculated_nps[step_idx], "\n",
                  "Calculated np uncertainty:", self.calculated_np_errors[step_idx],)

    def get_monte_carlo_step(self, pi, pf, pf_minus1, pcf_minus1, np_minus_1, np_minus_1_uncertainty):

        np_step_input_array = VaporSorptionSystem._create_np_input_step_array(
            t=self.temperature, pi=pi, pf=pf, pf_minus1=pf_minus1, pcf_minus1=pcf_minus1, np_minus1=np_minus_1)

        np_step_variance_array = VaporSorptionSystem._create_np_variance_array_for_experiment(
            input_step_array=np_step_input_array, np_minus1_variance=np_minus_1_uncertainty
        )

        np_uncertainty, np_inputs, np_mc_outputs = VaporSorptionSystem.get_monte_carlo_results(
            self.analytical_np_function, np_step_input_array, np_step_variance_array,
            iterations=VaporSorptionSystem.DEFAULT_ITERATIONS
        )
        np_analytical = self.analytical_np_function(*np_step_input_array)

        conc_step_input_array = VaporSorptionSystem._create_concentration_input_array(
            n_p=np_analytical, pen_mw=self.pen_mw, pol_dry_mass=self.polymer_mass, pol_dens=self.polymer_density)

        conc_step_variance_array = VaporSorptionSystem._create_concentration_variance_array_for_experiment(
            np_variance=np_uncertainty, pen_mw_variance=self.pen_mw_error,
            polymer_mass_variance=self.polymer_mass_error, polymer_density_variance=self.polymer_density_error)
        conc_uncertainty, conc_inputs, conc_mc_outputs = VaporSorptionSystem.get_monte_carlo_results(
            self.analytical_concentration_function, conc_step_input_array, conc_step_variance_array,
            iterations=VaporSorptionSystem.DEFAULT_ITERATIONS
        )
        conc_analytical = self.analytical_concentration_function(*conc_step_input_array)
        # print(VaporSorptionSystem.DEFAULT_ITERATIONS)
        # print("Conc_uncertainty:", conc_uncertainty)
        # print("Conc_value      :", conc_analytical)
        # print("Np_uncertainty  :", np_uncertainty)
        # print("Np_value        :", np_analytical)

        np_pressures = [np_input[2] for np_input in np_inputs]
        # plotstuff(np_pressures, conc_mc_outputs, conc_mc_outputs)

        return (np_analytical, np_uncertainty, np_mc_outputs), (conc_analytical, conc_uncertainty, conc_mc_outputs)

    @staticmethod
    def get_monte_carlo_results(function, input_array, variance_array, iterations):

        # concentration_function = VaporSorptionSystem._get_analytical_concentration_function()
        # np_function = VaporSorptionSystem._get_analytical_np_function(pressure_dependent_error=True)

        function_inputs, function_outputs = ErrorPropagation.monte_carlo_error(
            model=function,
            input_array=input_array,
            variance_array=variance_array,
            iterations=iterations,
            normal_distribution=VaporSorptionSystem.USE_NORMAL_DISTRIBUTION)

        np_uncertainty = ErrorPropagation.get_uncertainty_from_monte_carlo_output(function_outputs)
        return np_uncertainty, function_inputs, function_outputs

    # Function generation

    @staticmethod
    def _get_analytical_np_function(pressure_dependent_error=True):
        """
        Function signature: f(R, T, pi, pf, vc, vs, pf_minus1, pcf_minus1, np_minus1)
            R -> J/molk
            T -> K
            pressure -> pascal
            vol -> cm3
        """
        npfunc = _VaporSorptionEquations(pressure_dependent_error=pressure_dependent_error).get_n_p_function()
        optimalnpfunc = nb.njit(npfunc)
        return optimalnpfunc

    @staticmethod
    def _get_analytical_concentration_function():
        """
        Function signature: f(n_p, pen_mw, pol_dry_mass, pol_dens)
        """
        # right now, this equation isn't dependent on whether or not the error is pressure dependent
        concfunc = _VaporSorptionEquations(
            pressure_dependent_error=True).get_concentration_function()
        optimalconcfunc = nb.njit(concfunc)
        return optimalconcfunc

    # LaTeX generation
    @staticmethod
    def get_analytical_np_variance_latex(pressure_dependent_error=True):
        return _VaporSorptionEquations(pressure_dependent_error=pressure_dependent_error).n_p_variance()

    @staticmethod
    def get_analytical_concentration_variance_latex():
        # right now, this equation isn't dependent on whether or not the error is pressure dependent
        return _VaporSorptionEquations(pressure_dependent_error=True).concentration_variance()

    # Helper methods
    @staticmethod
    def _create_np_variance_array_for_experiment(input_step_array, np_minus1_variance):

        # signature: (R, T, pi, pf, vc, vs, pf_minus1, pcf_minus1, np_minus1)

        pi_measurement = input_step_array[2]
        pf_measurement = input_step_array[3]
        pf_minus1_measurement = input_step_array[6]
        pcf_minus1_measurement = input_step_array[7]

        temp_var = VaporSorptionSystem.TEMPERATURE_MEASUREMENT_VARIANCE  # K
        pres_var = VaporSorptionSystem.PRESSURE_MEASUREMENT_VARIANCE  # percentage
        vc_var = VaporSorptionSystem.VC_VARIANCE
        vs_var = VaporSorptionSystem.VS_VARIANCE
        # signature: (R, T, pi, pf, vc, vs, pf_minus1, pcf_minus1, np_minus1)
        variance_step = [0, temp_var, pi_measurement * pres_var, pf_measurement * pres_var,
                         vc_var, vs_var, pf_minus1_measurement * pres_var, pcf_minus1_measurement * pres_var,
                         np_minus1_variance]

        return variance_step

    @staticmethod
    def _create_concentration_variance_array_for_experiment(
            np_variance, pen_mw_variance, polymer_mass_variance, polymer_density_variance):
        # signature: (n_p, pen_mw, pol_dry_mass, pol_dens)
        return [np_variance, pen_mw_variance, polymer_mass_variance, polymer_density_variance]

    @staticmethod
    def _create_np_input_step_array(t, pi, pf, pf_minus1, pcf_minus1, np_minus1):
        # signature: (R, T, pi, pf, vc, vs, pf_minus1, pcf_minus1, np_minus1)
        return [VaporSorptionSystem.R, t, pi, pf, VaporSorptionSystem.VC, VaporSorptionSystem.VS,
                pf_minus1, pcf_minus1, np_minus1]

    @staticmethod
    def _create_concentration_input_array(n_p, pen_mw, pol_dry_mass, pol_dens):
        # signature: (n_p, pen_mw, pol_dry_mass, pol_dens)
        return [n_p, pen_mw, pol_dry_mass, pol_dens]

    @staticmethod
    def test():
        test_passed = True
        sigfigs = 10

        # check vapor sorption equations for consistency
        result = VaporSorptionSystem._get_analytical_concentration_function()(6.1437e-5, 32, 0.0234, 1.393)
        expected = 81.97560439205128
        if not round(expected, sigfigs) == round(result, sigfigs):
            test_passed = False
            print("Analytical concentration function failed numerical check. Expected:", expected, "Result:", result)
        result = VaporSorptionSystem._get_analytical_np_function(
            pressure_dependent_error=True)(
            8.314, 298.15, 5857.47381578947, 3801.41898538673, 29.4778313316133,
            7.61387679765219, 1860.51134470043, 1867.33029499151, 4.30292156079647e-5)
        expected = 6.143685091334635e-05
        if not round(expected, sigfigs) == round(result, sigfigs):
            test_passed = False
            print("Analytical moles sorbed function failed numerical check. Expected:", expected, "Result:", result)

        # check single step in the VaporSorptionSystem
        sorption_system = VaporSorptionSystem(
            polymer_density=1.393, polymer_density_error=0.001,
            polymer_mass=0.0234, polymer_mass_error=5e-5,
            penetrant_mw=32, pen_mw_error=0,
            temperature=298.15, temperature_error=1,
        )

        # Deprecated
        # print(sorption_system.get_monte_carlo_step(
        #     pi=2681.81466294, pf=643.2713876919, pf_minus1=0, pcf_minus1=0,
        #     np_minus_1=0, np_minus_1_uncertainty=0))

        # actually use the sorption system in the correct manner
        sorption_system.add_step_data(pi=2681.81466294, pf=643.27138769,  pcf=643.947039477)
        sorption_system.add_step_data(pi=3921.56620664, pf=1860.51134470, pcf=1867.33029499)
        sorption_system.add_step_data(pi=5857.47381579, pf=3801.41898539, pcf=3818.82878289)
        sorption_system.add_step_data(pi=7107.19325658, pf=5635.72490346, pcf=5643.36920230)
        sorption_system.add_step_data(pi=8091.36872470, pf=6991.55832237, pcf=7025.70789474)
        sorption_system.add_step_data(pi=8091.36872470, pf=6991.55832237, pcf=7025.70789474)
        sorption_system.add_step_data(pi=9556.43143593, pf=8417.08552632, pcf=None)
        print(sorption_system.get_monte_carlo_isotherm())
        # check the monte carlo single step calculation for errors / surface level behavior
        # np_input_step_1 = VaporSorptionSystem._create_np_input_step_array(
        #     t=298.15, pi=2681.81466294, pf=643.2713876919, pf_minus1=0, pcf_minus1=0, np_minus1=0)
        #
        # np_variance_step_1 = VaporSorptionSystem._create_np_variance_array_for_experiment(
        #     input_step_array=np_input_step_1, np_minus1_variance=0.00  # step one has no np_minus1, so no need to vary
        # )
        # conc_input_step_1 = VaporSorptionSystem._create_concentration_input_array(
        #     np=,
        #
        # )
        # conc_variance_step_1 = VaporSorptionSystem._create_concentration_variance_array_for_experiment(
        #
        # )
        # np_uncertainty_step_1, conc_uncertainty_step_1 = VaporSorptionSystem.get_monte_carlo_results(
        #     np_input_array=np_input_step_1,
        #     np_variance_array=np_variance_step_1,
        #     concentration_input_array=,
        #     concentration_variance_array=,
        #     iterations=1000000))

        return test_passed

# VaporSorptionSystem().concentration_variance()
# VaporSorptionSystem().n_p_variance()
# VaporSorptionSystem(pressure_dependent_error=False).n_p_variance()


if __name__ == "__main__":
    if VaporSorptionSystem.test():
        print("All tests passed.")

    print("Moles absorbed variance:", VaporSorptionSystem.get_analytical_np_variance_latex(
        pressure_dependent_error=True))
    print("Concentration variance :", VaporSorptionSystem.get_analytical_concentration_variance_latex())
