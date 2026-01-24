"""Module alias for orthogonal polynomial recurrences."""
try:
    from recursum._recursum import (
        legendre,
        chebyshevt,
        chebyshevu,
        hermite,
        hermiteh,
        hermitehe,
        laguerre,
        assoclegendre,
    )
    __all__ = [
        'legendre',
        'chebyshevt',
        'chebyshevu',
        'hermite',
        'hermiteh',
        'hermitehe',
        'laguerre',
        'assoclegendre',
    ]
except ImportError:
    pass
