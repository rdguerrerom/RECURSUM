"""Module alias for special function recurrences."""
try:
    from recursum._recursum import (
        jacobi,
        gegenbauer,
        assoclaguerre,
        airyai,
        airybi,
        besselj,
        bessely,
        sphericalbesselj,
        sphericalbessely,
        modifiedbesseli,
        modifiedbesselk,
        euler,
        bernoulli,
    )
    __all__ = [
        'jacobi',
        'gegenbauer',
        'assoclaguerre',
        'airyai',
        'airybi',
        'besselj',
        'bessely',
        'sphericalbesselj',
        'sphericalbessely',
        'modifiedbesseli',
        'modifiedbesselk',
        'euler',
        'bernoulli',
    ]
except ImportError:
    pass
