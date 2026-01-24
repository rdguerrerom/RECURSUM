"""Module alias for McMurchie-Davidson recurrences."""
try:
    from recursum._recursum import (
        hermitee,
    )
    __all__ = ['hermitee']
except ImportError:
    pass
