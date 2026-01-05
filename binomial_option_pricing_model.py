# TODO: Numba optimization
# There is a possibility to improve speed of following code by using numba library

# Import core Python libraries
import time
import math


def price_option(
        spot: int | float,
        strike: int | float,
        sigma: int | float,
        steps: int,
        T: int | float,
        rate: int | float,
        side: str | int = 'call',
        kind: str | int = 'european',
        relative: bool = False
) -> tuple[float, float, int]:

    '''
    This function calculates the option price using the binomial tree model.

    :param spot: current spot price of the underlying asset
    :param strike: strike price of the option
    :param sigma: volatility (standard deviation) of the underlying asset
    :param steps: number of steps in the binomial tree
    :param T: time to expiration of the option in years
    :param rate: current risk-free rate
    :param side: option side (call or 1, put or -1), default: call
    :param kind: option type (european or 0, american or 1), default: european
    :param relative: sets time measurement to be absolute (real time) or relative (for performance comparison)

    :return: option price, computation time, number of paths considered
    '''

    # Start measuring computation time
    if relative:
        start = time.perf_counter()
    else:
        start = time.time()

    # Validate input parameters
    if not isinstance(spot, (int, float)):
        raise TypeError(f"Expected int or float, but got {type(spot).__name__}")
    if spot <= 0:
        raise ValueError("Spot needs to be grater than 0")
    if not isinstance(strike, (int, float)):
        raise TypeError(f"Expected int or float, but got {type(strike).__name__}")
    if strike <= 0:
        raise ValueError("Strike needs to be grater than 0")
    if not isinstance(sigma, (int, float)):
        raise TypeError(f"Expected int or float, but got {type(sigma).__name__}")
    if sigma <= 0:
        raise ValueError("Standard deviation needs to be grater than 0")
    if not isinstance(steps, int):
        raise TypeError(f"Expected int, but got {type(steps).__name__}")
    if steps <= 0:
        raise ValueError("Steps need to be grater than 0")
    if not isinstance(T, (int, float)):
        raise TypeError(f"Expected int, float, but got {type(T).__name__}")
    if T <= 0:
        raise ValueError("Time to maturity needs to be greater than 0")
    if not isinstance(rate, (int, float)):
        raise TypeError(f"Expected int, float, but got {type(T).__name__}")
    if not isinstance(side, (str, int)):
        raise TypeError(f"Expected int, str, but got {type(T).__name__}")
    if side not in {'call', 'put', -1, 1}:
        raise ValueError("Unexpected value, try 'call' or 'put', or 1 (call) or -1 (put)")
    if not isinstance(kind, (str, int)):
        raise TypeError(f"Expected int, str, but got {type(T).__name__}")
    if kind not in {'european', 'american', 0, 1}:
        raise ValueError("Unexpected value, try 'european' or 'american', or 0 (european) or 1 (american)")

    # Set numerical flags for side and option type
    if type(side) is not int:
        side = 1 if side == 'call' else -1
    if type(kind) is not int:
        kind = int(kind != 'european')

    # Calculate constants

    # Time increment per step
    t = T / steps
    # Price diffusion coefficient
    diffusion = sigma * math.sqrt(t)

    # Up and down movement factors
    u = math.exp(diffusion)
    # u * d = exp(diffusion) * exp(-diffusion) = 1
    d = 1.0 / u

    # Using exponentiation is computationally expensive, so we use multiplication
    # to move from one terminal node to the next.
    # The direction depends on the option side. See README section 3.1.
    move = d * d if side == 1 else u * u

    # Calculate risk-neutral probabilities
    p = (math.exp(rate * t) - d) / (u - d)
    q = 1.0 - p

    # Calculate discount factor per step
    df = math.exp(-rate * t)

    # Calculate the number of paths considered
    considered_steps = 2**steps

    # Initialize payoff list
    payoffs = [0.0] * (steps + 1)

    # Calculate the index of the terminal node at which the option becomes worthless.
    # See README section 3.2.
    k = math.ceil(steps / 2 - (math.log(strike / spot if side == 1 else spot / strike)) / (2 * diffusion))
    # Cap the index at the number of steps (if k > steps, all nodes are in the money)
    k = min(steps + 1, k)

    # If k is less than or equal to zero, the option never reaches the strike price
    if k <= 0:
        return 0.0, time.time() - start, considered_steps

    # Calculate the most extreme terminal price based on option side
    curr_spot = spot * (u ** steps if side == 1 else d ** steps)

    # Calculate terminal payoffs
    for i in range(k):
        # Call: S - K, Put: K - S = (S - K) * (-1)
        val = side * (curr_spot - strike)
        # Assign payoff to the corresponding node
        # Using conditional checks is faster than max()
        payoffs[i] = val if val > 0 else 0.0
        # Move to the next terminal node price
        curr_spot *= move

    # This loop is only used for American puts.
    # See README section 3.3.
    if kind and side == -1:
        # American option pricing loop
        # Iterate backwards from the n-1 step to the root
        for i in range(steps - 1, -1, -1):

            # Start at the lowest node for this time step: S * d^i
            # Puts are deepest in-the-money at low prices
            node_spot = spot * (d ** i)

            for j in range(min(i + 1, k)):
                # Expected value
                expected_payoff = df * (p * payoffs[j + 1] + q * payoffs[j])

                # Intrinsic value
                intrinsic = strike - node_spot

                # Choose between early exercise and continuation
                payoffs[j] = intrinsic if intrinsic >= expected_payoff and intrinsic > 0 else (expected_payoff if expected_payoff > 0 else 0)

                # Move spot price up to the next node
                node_spot *= move

    # European options and American calls
    else:

        # To optimize performance, we use the binomial distribution
        # in logarithmic form to avoid numerical overflow.
        # If we try to use n!/((n-j)!j!), with large enough steps variable
        # we get OverflowError: integer division result too large for a float
        # See README section 3.4.

        log_p = math.log(p)
        log_q = math.log(q)
        log_fact_n = math.lgamma(steps + 1)
        expected_val = 0.0

        # Call option
        if side == 1:

            for j in range(k):
                # Compute binomial coefficient in log space
                log_comb = log_fact_n - math.lgamma(j + 1) - math.lgamma(steps - j + 1)
                # Compute the probability of reaching node j
                log_prob = log_comb + (steps - j) * log_p + j * log_q
                # Add the weighted payoff to the expected value
                expected_val += math.exp(log_prob) * payoffs[j]

        # Put option
        # Same logic as above, but with reversed probabilities
        else:
            for j in range(k):
                log_comb = log_fact_n - math.lgamma(j + 1) - math.lgamma(steps - j + 1)
                log_prob = log_comb + j * log_p + (steps - j) * log_q
                expected_val += math.exp(log_prob) * payoffs[j]

        # Discount expected value to present value
        payoffs[0] = math.exp(-rate * T) * expected_val


    if relative:
        return payoffs[0], time.perf_counter() - start, considered_steps
    else:
        return payoffs[0], time.time() - start, considered_steps


if __name__ == "__main__":
    # Import library to handle large number of paths
    import decimal

    # Set parameters
    spot = 10
    strike = 5
    sigma = 0.5
    steps = 10_000
    T = 1
    rate = 0.01
    side = "call"
    kind = "european"

    # Run the function
    option_price, computation_time, number_of_paths = price_option(spot, strike, sigma, steps, T, rate, side, kind)

    # Helper function to convert computation time into readable units
    def convert_time(time_to_convert: int | float) -> str:
        '''
        :param time_to_convert: time in seconds
        :return: formatted string with time in microseconds, milliseconds, or seconds
        '''

        if not isinstance(time_to_convert, (int, float)):
            raise TypeError(f"Expected int, float, but got {type(time_to_convert).__name__}")

        if time_to_convert < 1e-3:
            return f"{time_to_convert * 1e6:.2f} μs"
        elif time_to_convert < 1:
            return f"{time_to_convert * 1e3:.2f} ms"
        else:
            return f"{time_to_convert:.2f} s"

    # Display function inputs
    print("\nFunction inputs\n"+ 50*"=")
    print(f"Spot price: {spot:.2f}\n"
          f"Strike price: {strike:.2f}\n"
          f"Volatility: {sigma:.2f}\n"
          f"Time till maturity: {T:.2f}\n"
          f"Steps: {steps:,d}\n"
          f"RFR: {rate*100:.2f}%\n"
          f"Side: {side}\n"
          f"Type: {kind}")

    # Display function outputs
    print("\nFunction outputs\n" + 50*"=")
    print(f"Option price: {option_price:.4f}\n"
          f"Computation time: {convert_time(computation_time)}\n"
          f"Paths considered: {format(decimal.Decimal(number_of_paths), '.4e')}\n")