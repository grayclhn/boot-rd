I don’t know Stata well enough to try to use their notation... so, sorry.

I think we want to start with a few functions.

- `lhs_polynomial`: just estimate a given polynomial local to the threshold,
  but only on the left side of the threshold. We can mirror this function
  with `rhs_polynomial` for the other half.

- `estimate_lhs_bias`: estimates the bias of the pth order local
  polynomial using a local bootstrap, for just the left side of the
  discontinuity

- `bias_corrected_estimator`: estimates the effect of the discontinuity
  using the bootstrap bias correction

- `bootstrap_distribtution`: calculates the distribution of the
  bias_corrected_estimator using the bootstrap.

Let’s also just use the uniform kernel for now. If the approach
doesn’t work at all for the uniform kernel, using a better kernel
won’t save it. We’re also assuming the discontinuity is at x = 0.

## Bias correction

I assume that there's a library of functions like this next one, so
we probably shouldn't have to use it. But if we do:

```
lhs_polynomial = function(y, x, p, h)
	# y is a vector of N observations on the outcome of interest
	# x is a vector of N observations on the running variable
	# p is the order of the polynomial that the researcher wants to fit
	# h is the width of the uniform kernel that we’ll use to fit the polynomial

	L = indices of x such that -h <= x < 0 # ‘L’ stands for ‘local’
	estimate a polynomial of order p using y[L] and x[L] and store it as g
  yHat = g(x[L]) # fitted values
  epHat = y[L] - yHat
	return g, yHat, epHat
end function

rhs_polynomial(y, x, p, h) = lhs_polynomial(y, -x, p, h)
```

Then we actually want to use the bootstrap to estimate the one-sided bias.

```
estimate_lhs_bias = function(y, x, p, h1, h2, nboots)
	# y, x, and p are defined as in `lhs_polynomial`
	# h1 is the smaller bandwidth that is used to fit the pth order polynomial
  #    of interest
	# h2 is the larger bandwidth that’s used to fit the approximating higher
  #    order polynomial
	# nboots is the number of bootstrap simulations to run

	L = indices of x such that -h2 <= x < 0
	M = number of elements in L

	gHigher, yHat, epHat = lhs_polynomial(y, x, p+1, h2)
	boot_intercept = gHigher(0)

	allocate a numeric vector of length nboots and call it `discrepancy`
	for b in 1 to nboots:
		epBoot = M values randomly drawn from epHat with replacement
		yBoot = yHat + epBoot
		gBoot = lhs_polynomial(y, x, p, h1)
		discrepancy[b] = gBoot(0) - boot_intercept
	end loop

	return mean(discrepancy) # as the estimate of the bias
end function

estimate_rhs_bias(y, x, p, h1, h2, nboots) =
  estimate_lhs_bias(y, -x, p, h1, h2, nboots)
```

Then the first main (but simple) function

```
bias_corrected_estimator = function(y, x, p, h1, h2, nboots)
	# Same variable definitions as for `estimate_lhs_bias`

	g_lhs = lhs_polynomial(y, x, p, h1)
	g_rhs = rhs_polynomial(y, x, p, h1)

	return (g_rhs(0) - estimate_rhs_bias(y, x, p, h1, h2, nboots)) -
	  (g_lhs(0) - estimate_lhs_bias(y, x, p, h1, h2, nboots)
end function
```

## Coverage

We probably (definitely) want to studentize these estimates, but... if it
doesn't have approximately correct coverage for the unstudentized versions,
it's not going to be saved by studentization. So let's put it off.

```
bootstrap_distribution = function(y, x, p, h1, h2, nboots1, nboots2)

  L_lhs = indices of x such that -h2 <= x < 0
  L_rhs = indices of x such that   0 < x <= h2
  M_lhs = number of elements in L_lhs
  M_rhs = number of elements in L_rhs

  lhs_higher, yHat_lhs, epHat_lhs = lhs_polynomial(y, x, p+1, h2)
  rhs_higher, yHat_rhs, epHat_rhs = rhs_polynomial(y, x, p+1, h2)

  allocate a numeric vector of length nboots1 called "boot_estimates"
  for b in 1:nboots1
    epBoot_lhs = draw length(epHat_lhs) values from epHat_lhs with replacement
    epBoot_rhs = same, from epHat_rhs
    yBoot_lhs = epBoot_lhs + yHat_lhs
    yBoot_rhs = epBoot_rhs + yHat_rhs
    boot_estimates[b] = bias_corrected_estimator([yBoot_lhs, yBoot_rhs],
      [x[L_lhs], x[L_rhs]], p, h1, h2, nboots2)
  end loop
  return boot_estimates
end function
```

## Monte Carlo

Then we need to actually run the simulations. This I'm just going to sketch out.

```
for i in 1:nsims
  generate x and y data from DGP
  select bandwidths h1 and h2 maybe based on data # is there an easy non-optimal way?
  boot_estimates = bootstrap_distribution(y, x, 1, h1, h2, 19, 19) # increase once it runs!
  interval = quantile(boot_estimates, 0.025, 0.975)
  check and store whether interval contains true value
end for
calculate coverage in Monte Carlo.
```
