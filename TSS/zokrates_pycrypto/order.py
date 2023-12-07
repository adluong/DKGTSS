"""
This code is copied from https://github.com/ethereum/py_ecc/blob/master/py_ecc/bn128/bn128_curve.py
Author is Vitalik Buterin.
Unfortunately the field modulus is not generic in this implementation, hence we had to copy the file.
All changes from our side are denoted with #CHANGE.
"""

from __future__ import absolute_import

from typing import cast, List, Tuple, Sequence, Union


# The prime modulus of the field
# field_modulus = 21888242871839275222246405745257275088696311157297823662689037894645226208583
field_modulus = (
    21888242871839275222246405745257275088614511777268538073601725287587578984328
)
# CHANGE: Changing the modulus to the embedded curve

# See, it's prime!
# assert pow(2, field_modulus, field_modulus) == 2

# The modulus of the polynomial in this representation of FQO12
# FQO12_MODULUS_COEFFS = (82, 0, 0, 0, 0, 0, -18, 0, 0, 0, 0, 0)  # Implied + [1]
# FQO2_MODULUS_COEFFS = (1, 0)
# CHANGE: No need for extended  in this case

# Extended euclidean algorithm to find modular inverses for
# integers
def inv(a: int, n: int) -> int:
    if a == 0:
        return 0
    lm, hm = 1, 0
    num = a if isinstance(a, int) else a.n
    low, high = num % n, n
    while low > 1:
        r = high // low
        nm, new = hm - lm * r, high - low * r
        lm, low, hm, high = nm, new, lm, low
    return lm % n


IntOrFQO = Union[int, "FQO"]


# A class for field elements in FQO. Wrap a number in this class,
# and it becomes a field element.
class FQO(object):
    n = None  # type: int

    def __init__(self, val: IntOrFQO) -> None:
        if isinstance(val, FQO):
            self.n = val.n
        else:
            self.n = val % field_modulus
        assert isinstance(self.n, int)

    def __add__(self, other: IntOrFQO) -> "FQO":
        on = other.n if isinstance(other, FQO) else other
        return FQO((self.n + on) % field_modulus)

    def __mul__(self, other: IntOrFQO) -> "FQO":
        on = other.n if isinstance(other, FQO) else other
        return FQO((self.n * on) % field_modulus)

    def __rmul__(self, other: IntOrFQO) -> "FQO":
        return self * other

    def __radd__(self, other: IntOrFQO) -> "FQO":
        return self + other

    def __rsub__(self, other: IntOrFQO) -> "FQO":
        on = other.n if isinstance(other, FQO) else other
        return FQO((on - self.n) % field_modulus)

    def __sub__(self, other: IntOrFQO) -> "FQO":
        on = other.n if isinstance(other, FQO) else other
        return FQO((self.n - on) % field_modulus)

    def __div__(self, other: IntOrFQO) -> "FQO":
        on = other.n if isinstance(other, FQO) else other
        assert isinstance(on, int)
        return FQO(self.n * inv(on, field_modulus) % field_modulus)

    def __truediv__(self, other: IntOrFQO) -> "FQO":
        return self.__div__(other)

    def __rdiv__(self, other: IntOrFQO) -> "FQO":
        on = other.n if isinstance(other, FQO) else other
        assert isinstance(on, int), on
        return FQO(inv(self.n, field_modulus) * on % field_modulus)

    def __rtruediv__(self, other: IntOrFQO) -> "FQO":
        return self.__rdiv__(other)

    def __pow__(self, other: int) -> "FQO":
        if other == 0:
            return FQO(1)
        elif other == 1:
            return FQO(self.n)
        elif other % 2 == 0:
            return (self * self) ** (other // 2)
        else:
            return ((self * self) ** int(other // 2)) * self

    def __eq__(
        self, other: IntOrFQO
    ) -> bool:  # type:ignore # https://github.com/python/mypy/issues/2783 # noqa: E501
        if isinstance(other, FQO):
            return self.n == other.n
        else:
            return self.n == other

    def __ne__(
        self, other: IntOrFQO
    ) -> bool:  # type:ignore # https://github.com/python/mypy/issues/2783 # noqa: E501
        return not self == other

    def __neg__(self) -> "FQO":
        return FQO(-self.n)

    def __repr__(self) -> str:
        return repr(self.n)

    def __int__(self) -> int:
        return self.n

    @classmethod
    def one(cls) -> "FQO":
        return cls(1)

    @classmethod
    def zero(cls) -> "FQO":
        return cls(0)
