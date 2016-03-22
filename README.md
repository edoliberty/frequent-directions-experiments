# Matrix Sketching
This repo was created by [Edo Liberty](www.edoliberty.com) and [Mina Ghashami](http://www.cs.utah.edu/~ghashami/).
It builds all common streaming matrix sketching algroithms in Python.
It is developed for academic use only and for reproducability of the results in the following papers:


[Simple and Deterministic Matrix Sketches: Edo Liberty](http://www.cs.yale.edu/homes/el327/papers/simpleMatrixSketching.pdf)

[Relative Errors for Deterministic Low-Rank Matrix Approximations: Mina Ghashami, Jeff M. Phillips](http://www.cs.utah.edu/~ghashami/papers/relative_err_soda.pdf)

[Frequent Directions: Simple and Deterministic Matrix Sketching: Mina Ghashami, Edo Liberty, Jeff M. Phillips, David P. Woodruff](http://www.cs.utah.edu/~ghashami/papers/fd_journal.pdf)

[Efficient Frequent Directions Algorithm for Sparse Matrices: Mina Ghashami, Edo Liberty, Jeff M. Phillips](http://arxiv.org/abs/1602.00412)


 + Usage
If you are only using the library, you will noly need to the "python" folder.
It contains an exampleUsage.py file for your convenience.

 
 + Running tests and experiments 
Running tests requires using the -m flag which is standard in python unittesting. 
For example, to run the bruteForce sketcher test, go to the parent directory (outside frequentdirection/) and run
```
python -m frequentdirections.test.testBruteForce
```

 + Contributing
Please feel free to send me pull requests. The test package is minimal. 
So, if you make chages to the core classes. Please also include the tests to cover your changes. 