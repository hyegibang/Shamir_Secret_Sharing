from random import randint

from numpy.polynomial.polynomial import Polynomial as Poly
import numpy.polynomial.polynomial as polynomial


class SSS(object):
    """
    This class serves as an implementation of Shamir's Secret Sharing scheme,
    which provides methods for managing shared secrets
    """
    def __init__(self, secret, n, k, p):
        """
        S: secret
        n: total number of shares
        k: recovery threshold
        p: prime, where p > S and p > n
        """

        self.n = n
        self.k = k
        self.p = p

        mBytes = secret.encode("utf-8")
        Sval = int(mBytes.encode('hex'), 16)
        self.S = Sval

        production_coefs = [Sval]
        for coef in range(k-1):
        	production_coefs.append(randint(1,p-1))

        self.production_poly = Poly(production_coefs)

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
    def _lagrange_interpolate(x, x_s, y_s, p):
        '''
        Find the y-value for the given x, given n (x, y) points;
        k points will define a polynomial of up to kth order.
        '''
        k = len(x_s)
        assert k == len(set(x_s)), "points must be distinct"
        def PI(vals):  # upper-case PI -- product of inputs
            accum = 1
            for v in vals:
                accum *= v
            return accum
        nums = []  # avoid inexact division
        dens = []
        for i in range(k):
            others = list(x_s)
            cur = others.pop(i)
            nums.append(PI(x - o for o in others))
            dens.append(PI(cur - o for o in others))
        den = PI(dens)
        num = sum([SSS._divmod(nums[i] * den * y_s[i] % p, dens[i], p)
                   for i in range(k)])
        return (SSS._divmod(num, den, p) + p) % p

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
        Reconstructs a shared secret, given at least self.k of the proper shares
        """

        if len(shares) < self.k:
            raise Exception("Need more participants")

        x = [a for a, b in shares]
        y = [b for a, b in shares]

        return SSS._lagrange_interpolate(0, x, y, self.p)

    def decode_value(self,result):
        def to_bytes(n, length, endianess='big'):
            h = '%x' % n
            s = ('0'*(len(h) % 2) + h).zfill(length*2).decode('hex')
            return s if endianess == 'big' else s[::-1]
        mBytes2 = to_bytes(result,((result.bit_length() + 7) // 8))
        m2 = mBytes2.decode("utf-8")
        print("Reconstructed secret is:", m2)
        print(result == self.S)

def user_input_encode():
    secret = raw_input("Enter Your Secret: ")
    n, k = raw_input("Enter number of shares and threshold separate by space: ").split()
    n, k = int(n), int(k)

    return secret, n, k


if __name__ == "__main__":
    p = (2 ** 127) - 1

    secret, n, k = user_input_encode()

    print("Secret encoded is:", secret)
    sss = SSS(secret, n, k, p)
    y = sss.construct_shares()
    result = int(sss.reconstruct_secret(y[0:k]))

    sss.decode_value(result)
