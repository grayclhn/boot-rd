# Bootstrap Confidence Intervals for Sharp Regression Discontinuity Designs with the Uniform Kernel

This repository has the LaTeX and R source code for the paper,
"Bootstrap Confidence Intervals for Sharp Regression Discontinuity
Designs with the Uniform Kernel," by Otávio Bartalotti, Gray Calhoun,
and Yang He.

The datasets in `R/LM2007data` were originally downloaded from Doug
Miller's website,
<http://faculty.econ.ucdavis.edu/faculty/dlmiller/statafiles> and
were used for the paper,

-   Ludwig, Jens and Douglas L. Miller. 2007. “Does Head Start Improve
    Children’s Life Chances? Evidence from a Regression Discontinuity
    Design.” *Quarterly Journal of Economics* 122 (1):159– 208.

Researchers interested in using the data for their own work should get
it from there, since it is far more likely to be complete and current
than ours.

## Rerunning analysis and building the pdf

The R script `R/run_everything.R` runs all of the R code necessary for
the simulations and empirical exercise presented in the paper. It must
be run from the `R` subdirectory and it will call all of the other
scripts in that directory.

You can generate the pdf by calling `pdflatex` and `bibtex` on the
file `main_paper.tex` repeatedly as usual.

## Contact information

Please contact the authors if you have questions or find errors and bugs:

- Otávio Bartalotti: <bartalot@iastate.ed>
- Gray Calhoun: <gcalhoun@iastate.edu>
- Yang He: <yanghe@iastate.edu>

## License and copying

Copyright (c) 2016, Otávio Bartalotti, Gray Calhoun, and Yang He.

All of the R code in the `R` directory is licensed under the MIT
"Expat" License:

> Permission is hereby granted, free of charge, to any person
> obtaining a copy of this software and associated documentation files
> (the "Software"), to deal in the Software without restriction,
> including without limitation the rights to use, copy, modify, merge,
> publish, distribute, sublicense, and/or sell copies of the Software,
> and to permit persons to whom the Software is furnished to do so,
> subject to the following conditions:
>
> The above copyright notice and this permission notice shall be
> included in all copies or substantial portions of the Software.
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
> EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
> MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
> NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
> BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
> ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
> CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
> SOFTWARE.
