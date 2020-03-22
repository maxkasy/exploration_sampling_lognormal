library(tidyverse)

k=6
n=100

test_data= tibble(
    D=sample(k, n, replace=T),
    Y = exp(rnorm(n,D/10,.2))
)

write_csv(test_data, "test/test_data.csv")
