import sympy
from numba import njit, prange


class ErrorPropagation:
    def __init__(self):
        pass  # maybe later this can be something more object oriented, but thus far there is no need

    @staticmethod
    def create_analytical_variance_function(equation, variable_variance_dict):
        """
        :param equation: Sympy expression (not an equals() expression, rather, the "something" in f(x, y, ...)=something
        :param variable_variance_dict: dict of sympy symbols or expressions
                formatted as variables and their corresponding variances, i.e., {variable: variance, ...}
                All variables referenced in the equation MUST be present as the same object in
                    memory, or they will not be recognized and therefore differentiated/included.
        :type variable_variance_dict: dict
        :return: Sympy expression of the propagated error
        """
        variable_derivative_dict = {}
        for variable in variable_variance_dict.keys():
            variable_derivative_dict[variable] = sympy.diff(equation, variable)

        dummy = sympy.Symbol("dummy")  # I honestly, cannot remember the easier way to create an empty expression
        # this doesn't impact anything significantly, so lets just pretend its an initializer and forget I did this :)
        result = dummy

        for variable in variable_derivative_dict.keys():
            result = result + (variable_variance_dict[variable] ** 2) * (variable_derivative_dict[variable]**2)

        result = result - dummy  # an equally stupid way for me to do this
        result = sympy.sqrt(result)
        return result

    @staticmethod
    @njit(parallel=True)
    def create_monte_carlo_variance_function(
            model,
            inputs: list,
            variances: list,
            iterations=100000,
            normal_distribution=True):

        """
        :param model:
        :param inputs:
        :param variances:
        :param iterations: Number of iterations to use in the simulation
        :param normal_distribution: bool, whether to use a normal distribution or sample variances evenly
        :type inputs: tuple(float)
        :type variances: tuple(float)
        :type model: Numba JIT-compiled function
        :type inputs: [iterable]

        :return: (list, list) maximum and minimum outputs for the result of the function
        """
        def returned_function():

            for _ in prange(iterations):
                result = model(*inputs)


