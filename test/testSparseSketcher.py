import unittest
from ..sketch.sparseSketcher import SparseSketcher as Sketcher
from scipy.sparse import rand

class testSparseSketcher(unittest.TestCase):

  def test_running(self):
    n = 100
    d = 20
    ell = 5
    A = rand(n, d, density = 0.001, format = 'lil') 
    sketcher = Sketcher(d,ell)

    for v in A:
        sketcher.append(v)
    sketch = sketcher.get()
    self.assertEqual(sketch.shape,(ell,d))

if __name__ == '__main__':
    unittest.main()
