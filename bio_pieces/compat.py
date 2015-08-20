try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from future.builtins import map

try:
    import unittest2
except ImportError:
    import unittest


