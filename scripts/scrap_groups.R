# Scrap from 08_dunlopgr_groups.Rmd


**All data**
  
  Fitting this as a reduced model to test that more complex model that differentiates between treatments fits better.


```{r}
# Fit model
# Use multiple starting points to avoid finding a local but not global maxima
all.mixfit1 <- normalmixEM(log(growth.asv$k), lambda = .5, mu = c(-2.5, -1), sigma = 0.3)
all.mixfit2 <- normalmixEM(log(growth.asv$k), lambda = .5, mu = c(-2, -0.5), sigma = 0.3)
all.mixfit3 <- normalmixEM(log(growth.asv$k), lambda = .5, mu = c(-2.5, -1.5), sigma = 0.3)

# Results
summary(all.mixfit1)
summary(all.mixfit2)
summary(all.mixfit3)
# mu refers to distribution means
# sigma refers to standard deviation
# lambda refers to mixing weights (ie what proportion of data falls under distribution 1 and 2)

# Change sigma starting point
all.mixfit4 <- normalmixEM(log(growth.asv$k), lambda = .5, mu = c(-2.5, -1), sigma = 0.5)
all.mixfit5 <- normalmixEM(log(growth.asv$k), lambda = .5, mu = c(-2.5, -1), sigma = 0.1)

# Results
summary(all.mixfit1)
summary(all.mixfit4)
summary(all.mixfit5)
```

```{r}
# Visualize
data.frame(x = all.mixfit1$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 0.25, colour = "black", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(all.mixfit1$mu[1], all.mixfit1$sigma[1], lam = all.mixfit1$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(all.mixfit1$mu[2], all.mixfit1$sigma[2], lam = all.mixfit1$lambda[2]),
                colour = "blue", lwd = 1.5) +
  ylab("Density") +
  theme_test()

data.frame(x = all.mixfit1$x) %>%
  ggplot() +
  geom_density(aes(x), binwidth = 0.25, colour = "black", 
               fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(all.mixfit1$mu[1], all.mixfit1$sigma[1], lam = all.mixfit1$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(all.mixfit1$mu[2], all.mixfit1$sigma[2], lam = all.mixfit1$lambda[2]),
                colour = "blue", lwd = 1.5) +
  ylab("Density") +
  theme_test()
```


### Slow vs fast

NOTE: it seems like bootstrapping the CI might be the better way to establish this see Stahl 2012 WIRE Comp Stat an citation therein.



The best way I could find to cluster the data was to bisect the distributions in each treatment. This seems to be how other popular gaussian mixture model packages work as well, such as Mclust. This means I will be more liklely to detect false positive results if I try to perform hypothesis tests with clustered data from the same treatments (because they are exclusionary), so I will avoid this.

I wonder if there is a better way to cluster the data, that allows overlap between clusters?
  
  A work around is to use the mu, sigma, and lambda parameters estimated by the model to calculate a t-test statistic for hypothesis testing. I'll try this. Based on the model parmaters, I know that the mixtools function estimated distirbutions with equal variances, so I can use a regular Student's t-test.

I want to use this as evidence that two groups exist within the overall community. If I reject the null hypothesis, the distributions estimated by mix-tools cannot be meaningfully differentiated from one-another, and analysis based on separating the data by groups is less justifiable.

Student's t-test forumla:

t = (mu1 - mu2)/(sp*sqrt(sigma1^2/n + sigma2^2/n))

sp = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2))/(n1+n2-2))

df =  n1 + n2 – 2

List of critical values: https://www.gradecalculator.tech/t-table/


**Succesional water control**

```{r}
# Variables for formula
mu1 <- S17y.mixfit1$mu[1]
mu2 <- S17y.mixfit1$mu[2]
sigma1 <- S17y.mixfit1$sigma[1]
sigma2 <- S17y.mixfit1$sigma[2]
n1 <- round(S17y.mixfit1$lambda[1]*nrow(S17y), digits=0) # round to nearest integer
n2 <- round(S17y.mixfit1$lambda[2]*nrow(S17y), digits=0)

# Sp
sp = sqrt(((n1-1)*sigma1^2 + (n2-1)*sigma2^2)/(n1+n2-2))

# T-statistic
t = (mu1 - mu2)/(sp*sqrt(1/n1 + 1/n2))
t

# Calculate degrees of freedom
df <-   n1 + n2 - 2
df

# Critical value at df
crit <- 1.962
isTRUE(abs(t) > crit)
```

For successional water control, reject null hypothesis


**Succesional water control**

```{r}
# Treatment model and data
treatment <- S17n.mixfit1
treatment.data <- S17n

# Variables for formula
mu1 <- treatment$mu[1]
mu2 <- treatment$mu[2]
sigma1 <- treatment$sigma[1]
sigma2 <- treatment$sigma[2]
n1 <- floor(treatment$lambda[1]*nrow(treatment.data)) # round down to nearest integer
n2 <- floor(treatment$lambda[2]*nrow(treatment.data)) # round down to nearest integer

# Sp
sp = sqrt(((n1-1)*sigma1^2 + (n2-1)*sigma2^2)/(n1+n2-2))

# T-statistic
t = (mu1 - mu2)/(sp*sqrt(1/n1 + 1/n2))
t

# Calculate degrees of freedom
df <-   n1 + n2 - 2
df

# Critical value at df
crit <- 1.98
isTRUE(abs(t) > crit)
```

Reject null hypothesis.


**Cropped water control**

```{r}
# Treatment model and data
treatment <- C3n.mixfit1
treatment.data <- C3n

# Variables for formula
mu1 <- treatment$mu[1]
mu2 <- treatment$mu[2]
sigma1 <- treatment$sigma[1]
sigma2 <- treatment$sigma[2]
n1 <- floor(treatment$lambda[1]*nrow(treatment.data)) # round down to nearest integer
n2 <- floor(treatment$lambda[2]*nrow(treatment.data)) # round down to nearest integer

# Sp
sp = sqrt(((n1-1)*sigma1^2 + (n2-1)*sigma2^2)/(n1+n2-2))

# T-statistic
t = (mu1 - mu2)/(sp*sqrt(1/n1 + 1/n2))
t

# Calculate degrees of freedom
df <-   n1 + n2 - 2
df

# Critical value at df
crit <- 1.97
isTRUE(abs(t) > crit)
```

Reject null hypothesis.


**Successional C-amended**

```{r}
# Treatment model and data
treatment <- S17y.mixfit1
treatment.data <- S17y

# Variables for formula
mu1 <- treatment$mu[1]
mu2 <- treatment$mu[2]
sigma1 <- treatment$sigma[1]
sigma2 <- treatment$sigma[2]
n1 <- floor(treatment$lambda[1]*nrow(treatment.data)) # round down to nearest integer
n2 <- floor(treatment$lambda[2]*nrow(treatment.data)) # round down to nearest integer

# Sp
sp = sqrt(((n1-1)*sigma1^2 + (n2-1)*sigma2^2)/(n1+n2-2))

# T-statistic
t = (mu1 - mu2)/(sp*sqrt(1/n1 + 1/n2))
t

# Calculate degrees of freedom
df <-   n1 + n2 - 2
df

# Critical value at df
crit <- 1.98
isTRUE(abs(t) > crit)
```

Reject null hypothesis.


**Cropped C-amended**

```{r}
# Treatment model and data
treatment <- C3y.mixfit1
treatment.data <- C3y

# Variables for formula
mu1 <- treatment$mu[1]
mu2 <- treatment$mu[2]
sigma1 <- treatment$sigma[1]
sigma2 <- treatment$sigma[2]
n1 <- floor(treatment$lambda[1]*nrow(treatment.data)) # round down to nearest integer
n2 <- floor(treatment$lambda[2]*nrow(treatment.data)) # round down to nearest integer

# Sp
sp = sqrt(((n1-1)*sigma1^2 + (n2-1)*sigma2^2)/(n1+n2-2))

# T-statistic
t = (mu1 - mu2)/(sp*sqrt(1/n1 + 1/n2))
t

# Calculate degrees of freedom
df <-   n1 + n2 - 2
df

# Critical value at df
crit <- 1.97
isTRUE(abs(t) > crit)
```

Reject null hypothesis.


### 16S copy number

Hypothesis: slow group taxa will have lower 16S copy number than fast group taxa. 

Based on growth rate vs 16S copy umber correlations from before, it is more likely that "slow" taxa will vary widely in 16S copy number and fast taxa will not.

Incorporate groupings back into main dataframe:

BROKEN BELOW

```{r}
# Merge classifications with other variables
growth.paprica.asv.S17n <- growth.paprica.asv %>%
  filter(Soil=="S17" & Amendment=="N") %>%
  mutate(log_k = log(k)) %>%
  left_join(S17n.mixfit.post, by=c("log_k"="x", "Soil", "Amendment")) %>%
  select(everything(), -comp.1, -comp.2, -log_k) %>%
  unique() # not sure why but the merge created duplicate rows

growth.paprica.asv.C3n <- growth.paprica.asv %>%
  filter(Soil=="C3" & Amendment=="N") %>%
  mutate(log_k = log(k)) %>%
  left_join(C3n.mixfit.post, by=c("log_k"="x", "Soil", "Amendment")) %>%
  select(everything(), -comp.1, -comp.2, -log_k) %>%
  unique() # not sure why but the merge created duplicate rows

growth.paprica.asv.S17y <- growth.paprica.asv %>%
  filter(Soil=="S17" & Amendment=="Y") %>%
  mutate(log_k = log(k)) %>%
  left_join(S17y.mixfit.post, by=c("log_k"="x", "Soil", "Amendment")) %>%
  select(everything(), -comp.1, -comp.2, -log_k) %>%
  unique() # not sure why but the merge created duplicate rows

growth.paprica.asv.C3y <- growth.paprica.asv %>%
  filter(Soil=="C3" & Amendment=="Y") %>%
  mutate(log_k = log(k)) %>%
  left_join(C3y.mixfit.post, by=c("log_k"="x", "Soil", "Amendment")) %>%
  select(everything(), -comp.1, -comp.2, -log_k) %>%
  unique() # not sure why but the merge created duplicate rows

# Merge all
growth.paprica.asv.label <- bind_rows(growth.paprica.asv.S17n, growth.paprica.asv.C3n, growth.paprica.asv.S17y, growth.paprica.asv.C3y)
```

```{r}
growth.paprica.asv.label %>%
  ggplot(aes(x=Soil, y=log(n16S), color=label)) +
  geom_boxplot() +
  theme_test()
```

Not even going to bother testing that.


### Change in abundance

Hypothesis: Fast taxa experience greater change in abundance (boom-bust strategy)

```{r}
growth.paprica.asv.label %>%
  mutate(change_abund = end_abund - start_abund) %>%
  ggplot(aes(x=Soil, y=log(change_abund), color=label)) +
  geom_boxplot() +
  theme_test()
```

### Start day

Slow and fast groups might tend to experience lag more or less.

Hypothesis: Fast taxa experience less lag.

```{r}
growth.paprica.asv.label %>%
  ggplot(aes(x=Soil, y=log(start_day), color=label)) +
  geom_boxplot() +
  theme_test()
```

Seems like "groups" don't separate well based on other growth variables.