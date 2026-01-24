"""Module alias for combinatorics recurrences."""
try:
    from recursum._recursum import (
        binomial,
        fibonacci,
    )
    __all__ = ['binomial', 'fibonacci']
except ImportError:
    pass
