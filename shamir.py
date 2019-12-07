from random import randint

from numpy.polynomial.polynomial import Polynomial as Poly
import numpy.polynomial.polynomial as polynomial


class SSS(object):
    """
    This class serves as an implementation of Shamir's Secret Sharing scheme,
    which provides methods for managing shared secrets
    """

    @staticmethod
    def lagrangePolynomialBasis(j, x, k):
        """
        Create a Lagrange basis polynomial
        j = current index of basis
        x = array of x values of the points [x1, x2, x3, ..., xk]
        k = threshold number of shares
        """

        polys = [
            Poly([-1 * x[m], 1]) / Poly([x[j] - x[m]])
            for m in range(k) if m != j
        ]

        return reduce(lambda acc, p: acc * p, polys, 1)
    @staticmethod
    def lagrangePolynomial(x, y, k):
        """
        Create a linear combination of Lagrange basis polynomials
        """

        return sum([y[j] * SSS.lagrangePolynomialBasis(j, x, k) for j in range(k)])

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
        print(self.production_poly)

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

        return SSS.lagrangePolynomial(x, y, self.k).coef[0] % self.p

def to_bytes(n, length, endianess='big'):
    h = '%x' % n
    s = ('0'*(len(h) % 2) + h).zfill(length*2).decode('hex')
    return s if endianess == 'big' else s[::-1]

if __name__ == "__main__":
    p = (2 ** 31) - 1

    m = raw_input("Enter Your Secret: ")
    n, k = raw_input("Enter number of shares and threshold separate by space: ").split()
    n, k = int(n), int(k)
    mBytes = m.encode("utf-8")
    mInt = int(mBytes.encode('hex'), 16)

    # mBytes2 = to_bytes(result,((result.bit_length() + 7) // 8))
    # m2 = mBytes2.decode("utf-8")
    # print(mInt == m2)

    S = mInt

    print("Secret encoded is:", S)
    sss = SSS(S, n, k, p)
    y = sss.construct_shares()
    print((y[0:k]))
    result = int(sss.reconstruct_secret(y[0:k]))
    print(result)
    print(result == S)
