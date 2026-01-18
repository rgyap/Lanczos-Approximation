# Lanczos-Approximation (Java version)

A command-line Java program that calculates the coefficients of Lanczos' approximation of the Gamma function.

## How to use?

Compile the `Main.java` file then execute the output on the command line. It will then ask for two numbers:
- $k$, any positive integer.
- $g$, any nonnegative real number. 

Next, it will ask how many digits would you like to see for each number that will be outputted (with the help of `java.math.BigDecimal.setScale()`).

It will then output $k+1$ numbers, say $a'_1, ..., a'_k$ (in the same order of the output), which are computed depending on the values of $k$ and $g$. With that, for any $z\in\mathbb{C}$ such that $\Re(z)>-0.5$, we have the following approximation:

$$\Gamma(z+1) \approx \sqrt{2\pi} {\left(z+g+\tfrac{1}{2}\right)}^{z+\tfrac{1}{2}} e^{-\left(z+g+\tfrac{1}{2}\right)}\left(a'_0 + \frac{a'_1}{z+1} + {...} + \frac{a'_k}{z+k} \right)$$

For more technical details about how this works, please refer to the PDF file, which should also include a Python implementation.

