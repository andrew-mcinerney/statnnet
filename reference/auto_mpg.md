# Auto MPG data

A cleaned version of the UCI Auto MPG data. The six observations with
missing horsepower are omitted, and the car-name identifier is not
retained. Numeric variables remain in their original units. Cylinder
count and origin are factors so that term-level model selection and
grouped summaries retain their original interpretation.

## Usage

``` r
auto_mpg
```

## Format

A data frame with 392 rows and 8 variables:

- mpg:

  Fuel economy in miles per gallon.

- cylinders:

  Number of cylinders, stored as a factor.

- displacement:

  Engine displacement in cubic inches.

- horsepower:

  Engine horsepower.

- weight:

  Vehicle weight in pounds.

- acceleration:

  Time to accelerate from 0 to 60 mph, in seconds.

- model_year:

  Model year, represented by the final two digits.

- origin:

  Vehicle origin: USA, Europe, or Japan.

## Source

Quinlan, R. (1993). *Auto MPG* data set. UCI Machine Learning
Repository. [doi:10.24432/C5859H](https://doi.org/10.24432/C5859H) .
Distributed under the Creative Commons Attribution 4.0 International
licence.

## References

Quinlan, R. (1993). Auto MPG. UCI Machine Learning Repository.
[doi:10.24432/C5859H](https://doi.org/10.24432/C5859H) .
