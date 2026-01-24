"""Module alias for quantum chemistry recurrences."""
try:
    from recursum._recursum import (
        stoauxb,
        boys,
        gaunt,
    )
    __all__ = ['stoauxb', 'boys', 'gaunt']
except ImportError:
    pass
