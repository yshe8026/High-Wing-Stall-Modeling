import numpy as np

# Define the function yy = f(xx)
def f(xx):
    return (-94738.66870330359 * xx**4 +
            314836.3898499948 * xx**3 +
            -392340.3248280777 * xx**2 +
            217290.48582230246 * xx +
            -45125.05373275733)

# Define the derivative of the function f'(xx)
def f_prime(xx):
    return (-4 * 94738.66870330359 * xx**3 +
            3 * 314836.3898499948 * xx**2 +
            -2 * 392340.3248280777 * xx +
            217290.48582230246)

# Implementing the Newton-Raphson method
def newton_raphson(xx_initial, tolerance=1e-7, max_iterations=1000):
    xx = xx_initial
    for iteration in range(max_iterations):
        f_value = f(xx)
        f_prime_value = f_prime(xx)
        if abs(f_prime_value) < 1e-10:  # Avoid division by zero
            print("Derivative is too small.")
            break
        xx_next = xx - f_value / f_prime_value
        if abs(xx_next - xx) < tolerance:  # Check for convergence
            return xx_next, iteration
        xx = xx_next
    return xx, max_iterations

# Initial guess for xx
xx_initial = 0.89

# Finding the root
xx_root, iterations = newton_raphson(xx_initial)

# Display the result
print(xx_root)
print(iterations)
