"""Module alias for Rys quadrature recurrences."""
try:
    from recursum._recursum import (
        rys2d,
        ryshrr,
        rysvrrfull,
        ryspoly,
    )
    __all__ = ['rys2d', 'ryshrr', 'rysvrrfull', 'ryspoly']
except ImportError:
    pass
