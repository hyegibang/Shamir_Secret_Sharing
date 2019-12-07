from random import randint

from numpy.polynomial.polynomial import Polynomial as Poly
import numpy.polynomial.polynomial as polynomial


class SSS(object):
    """
    This class serves as an implementation of Shamir's Secret Sharing scheme,
    which provides methods for managing shared secrets
    """
    @staticmethod
    def _extended_gcd(a, b):
        '''
        Division in integers modulus p means finding the inverse of the
        denominator modulo p and then multiplying the numerator by this
        inverse (Note: inverse of A is B such that A*B % p == 1) this can
        be computed via extended Euclidean algorithm
        http://en.wikipedia.org/wiki/Modular_multiplicative_inverse#Computation
        '''
        x = 0
        last_x = 1
        y = 1
        last_y = 0
        while b != 0:
            quot = a // b
            a, b = b, a%b
            x, last_x = last_x - quot * x, x
            y, last_y = last_y - quot * y, y
        return last_x, last_y

    @staticmethod
    def _divmod(num, den, p):
        '''Compute num / den modulo prime p

        To explain what this means, the return value will be such that
        the following is true: den * _divmod(num, den, p) % p == num
        '''
        inv, _ = SSS._extended_gcd(den, p)
        return num * inv

    @staticmethod
    def _lagrange_interpolate(x_s, y_s, p):
        '''
        given n (x, y) points; k points will define a polynomial of up to kth order.
        and return the y-intercept of the polynomial
        '''

        # number of shares neededis how many points are inputted
        k = len(x_s)

        # Asserk that all of the points are distinct
        assert k == len(set(x_s)), "points must be distinct"

        # Function to calculate product of a series of numbers
        def PI(vals):
            accum = 1
            for v in vals:
                accum *= v
            return accum

        nums = []
        dens = []
        for i in range(k):
            others = list(x_s) # get all x values
            cur = others.pop(i) # pop off the one for the numerator
            # Manually create common denominator
            nums.append(PI(o for o in others)) # numerator is mult. of all but popped
            dens.append(PI(o - cur for o in others)) # mult. all denominators
        den = PI(dens) # mult. all denominators to get GCD
        num = sum([SSS._divmod(nums[i] * den * y_s[i] % p, dens[i], p)
                   for i in range(k)]) # multiply by common dinominator and weight by y before dividing each numerator by its corresponding denominator
        return (SSS._divmod(num, den, p) + p) % p # return lin

    def __init__(self, S, n, k, p):
        """
        S: secret
        n: total number of shares
        k: recovery threshold
        p: prime, where p > S and p > n
        """

        self.S = S
        self.n = n
        self.k = k
        self.p = p

        #production_coefs = [1234, 166, 94]
        production_coefs = [S]
        for coef in range(k-1):
        	production_coefs.append(randint(1,p-1))

        self.production_poly = Poly(production_coefs)

    def construct_shares(self):
        """
        Used to generate shares in a production environment, based on a
        known number of total shares
        """

        return [
            (x, polynomial.polyval(x, self.production_poly.coef) % self.p)
            for x in range(1, self.n + 1)
        ]

    def reconstruct_secret(self, shares):
        """
        Reconstructs a shared secret, given at least k of the proper shares
        """

        if len(shares) < self.k:
            raise Exception("Need more participants")

        x = [a for a, b in shares]
        y = [b for a, b in shares]

        x_s = x[0:k]
        y_s = y[0:k]

        return SSS._lagrange_interpolate(x_s, y_s, self.p)


def to_bytes(n, length, endianess='big'):
    h = '%x' % n
    s = ('0'*(len(h) % 2) + h).zfill(length*2).decode('hex')
    return s if endianess == 'big' else s[::-1]

if __name__ == "__main__":
    p = (2 ** 127) - 1

    m = raw_input("Enter Your Secret: ")
    n, k = raw_input("Enter number of shares and threshold separate by space: ").split()
    n, k = int(n), int(k)
    mBytes = m.encode("utf-8")
    mInt = int(mBytes.encode('hex'), 16)

    S = mInt

    print("Secret encoded is:", S)
    sss = SSS(S, n, k, p)
    y = sss.construct_shares()
    print("shares to reconstruct:", (y[0:k]))
    result = int(sss.reconstruct_secret(y[0:k]))
    print(result)


    mBytes2 = to_bytes(result,((result.bit_length() + 7) // 8))
    m2 = mBytes2.decode("utf-8")
    print(m2)
    print(result == S)
