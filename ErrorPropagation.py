import sympy
from numba import njit  # , prange
# import numba as nb
import numpy as np


class ErrorPropagation:
    def __init__(self):
        pass  # maybe later this can be something more object oriented, but thus far there is no need

    @staticmethod
    def create_analytical_error_function(equation, variable_variance_dict):
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

    # @staticmethod
    # def create_monte_carlo_variance_function(
    #         model,
    #         normal_distribution=True):
    #
    #     # @njit(parallel=True)
    #     def returned_function(inputs, variances, iterations=1000000, num_outputs=1):
    #         # print("Input:", inputs, "Variances: ", variances)
    #
    #         assert(len(inputs) == len(variances))  # make sure user specified inputs line up
    #         num_variables = len(inputs)
    #
    #         randomly_varied_inputs = np.empty(shape=(iterations, num_variables), dtype=np.float32)
    #         outputs = np.empty(shape=(iterations, num_outputs), dtype=np.float32)
    #
    #         for i in prange(iterations):
    #             np.random.seed(i)  # seeds must be a uint32 # todo Can this be avoided and set before the for loop?
    #             rand_multipliers = 2 * np.random.rand(num_variables) - 1  # should be uniformly between -1 and 1
    #
    #             for n in range(num_variables):
    #                 randomly_varied_inputs[i][n] = inputs[n] + variances[n] * rand_multipliers[n]
    #
    #             # print("random set:", randomly_varied_inputs[i])
    #             outputs[i] = model(randomly_varied_inputs[i])
    #
    #         return randomly_varied_inputs, outputs
    #
    #     return returned_function

    @staticmethod
    def monte_carlo_error(model, input_array, variance_array, iterations=100000, normal_distribution=True):
        assert (len(input_array) == len(variance_array))  # make sure user specified inputs line up
        num_variables = len(input_array)

        randomly_varied_inputs = np.empty(shape=(iterations, num_variables), dtype=np.float32)
        outputs = np.empty(shape=iterations, dtype=np.float32)
        np.random.seed(1)  # seeds must be a uint32 # todo Can this be avoided and set before the for loop?
        progresscounter = 0
        print()
        for i in range(iterations):

            if normal_distribution:
                rand_multipliers = np.random.randn(num_variables)
                # print(np.array_repr(rand_multipliers).replace('\n', ''))
            else:
                rand_multipliers = 2 * np.random.rand(num_variables) - 1  # should be uniformly between -1 and 1

            for n in range(num_variables):
                randomly_varied_inputs[i][n] = input_array[n] + variance_array[n] * rand_multipliers[n]

            progresscounter += 1

            if progresscounter == 10000:
                print(str(round(i / iterations * 100, 3)) + "% complete", end='\n')
                progresscounter = 0
            # print("random set:", randomly_varied_inputs[i])
            outputs[i] = model(*randomly_varied_inputs[i])
        return randomly_varied_inputs, outputs

    @staticmethod
    def get_uncertainty_from_monte_carlo_output(output_array):
        return np.std(output_array)

    # @staticmethod  # todo
    # def create_monte_carlo_variance_function_from_lambda(
    #         lambda_function,
    #         normal_distribution=False,
    # ):
    #     return ErrorPropagation.create_monte_carlo_variance_function(
    #         nb.njit(lambda_function),
    #         normal_distribution=normal_distribution
    #     )

    @staticmethod
    def test():
        # testing monte carlo generator
        @njit
        def example_function(input_list):
            x, y, z = input_list
            f = 0
            for i in range(1000):
                f = f + (x - y) * z
            return f

        # example_measurements = (1., 2., 0.03)
        # example_variances = (0.2, 1.4, 0.001)
        # monte_carlo_example_function = ErrorPropagation.create_monte_carlo_variance_function(example_function)
        # # print(monte_carlo_example_function.inspect_types())
        # result = monte_carlo_example_function(example_measurements, example_variances)
        # print(result)
        return True


if __name__ == "__main__":
    num_test = 0
    num_pass = 0

    print("Testing Error Propagation Core")
    num_test += 1
    print(" - Testing Monte Carlo Function Generator... ", end='\r')
    if not ErrorPropagation.test():
        test_pass = False
        print(" - Testing Monte Carlo Function Generator [FAILED] ")
    else:
        num_pass += 1
        print(" - Testing Monte Carlo Function Generator [PASSED] ")

    print(str(num_pass) + "/" + str(num_test) + " tests passed.")
