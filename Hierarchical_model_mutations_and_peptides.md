Hierarchical\_model\_mutations\_and\_peptides.Rmd
================
Jacqueline Buros Novik
7/5/2017

Summarize data
--------------

In this analysis we will be looking at how the number of `mutations` & `peptides` (predicted neoantigens) varies with the sample type (solid / ascites) and timing of acquisition (relapse / primary and treated / untreated).

``` r
ggplot(md, aes(x = mutations_per_mb, fill = specific_treatment)) + 
  facet_wrap(~tissue_type) +
  geom_histogram(position = 'dodge') +
  theme_minimal()
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](Hierarchical_model_mutations_and_peptides_files/figure-markdown_github-ascii_identifiers/plot-data-1.png)

How many observations do we have for each of these categories?

``` r
md %>%
  group_by(tissue_type, treatment, timepoint) %>%
  tally() %>%
  tidyr::spread(key = tissue_type, value = n, fill = 0)
```

    ## # A tibble: 3 x 4
    ## # Groups:   treatment [2]
    ##         treatment  timepoint ascites solid
    ## *           <chr>      <chr>   <dbl> <dbl>
    ## 1   chemo treated    primary       0     5
    ## 2   chemo treated recurrence      24     6
    ## 3 treatment naive    primary       4    75

Strikes me that the "recurrent" timepoint is problematic in this analysis, since we don't have any untreated/recurent samples & so cannot separate effect of recurrence from that of treatment.

Instead, we can look at `solid` samples only, among those collected at the primary timepoint comparing treated to untreated samples. My guess is, many of these might be the paired samples.

We can then include the treated / recurrence using only solid samples, although my guess is in this case we will see a higher rate of mutations among recurrence samples than among primary/treated.

Only after these effects are well established should we turn to the ascites samples, to see if the difference between untreated/primary & treated/relapse is consisent with that seen in solid samples.

Restricting to primary, solid samples
-------------------------------------

We now have samples, are treatment-naive.

Let's review how these metrics are distributed in this subset of our samples. First we note that the maximum number of samples per donor in this subset of our data is , meaning we have no duplicate samples.

Here, looking at metrics among all primary, solid samples irrespective of treatment:

``` r
ggplot(md_primary_solid %>%
         tidyr::gather(value = 'value', key = 'variable', mutations, mutations_per_mb, peptides), aes(x = value)) + 
  geom_density() + 
  theme_minimal() + 
  facet_wrap(~variable, scale = 'free')
```

![](Hierarchical_model_mutations_and_peptides_files/figure-markdown_github-ascii_identifiers/psolid-review-dist-1.png)

These numbers are not exactly normally-distributed.

Perhaps using a log-transformed value?

``` r
ggplot(md_primary_solid %>%
         tidyr::gather(value = 'value', key = 'variable', mutations, mutations_per_mb, peptides),
       aes(x = log1p(value))) + 
  geom_density() + 
  theme_minimal() + 
  facet_wrap(~variable, scale = 'free')
```

![](Hierarchical_model_mutations_and_peptides_files/figure-markdown_github-ascii_identifiers/psolid-review-dist-log-1.png)

``` r
ggplot(md_primary_solid %>%
         tidyr::gather(value = 'value', key = 'variable', mutations, mutations_per_mb, peptides),
       aes(x = log1p(value), fill = treatment)) + 
  geom_density(alpha = 0.4) + 
  theme_minimal() + 
  facet_wrap(~variable, scale = 'free')
```

![](Hierarchical_model_mutations_and_peptides_files/figure-markdown_github-ascii_identifiers/psolid-review-dist-log-by-trt-1.png)

What is noticeable here is that, given the small number of treated samples, it is very hard to tell graphically whether there is any difference in the two distributions.

Let's try fitting a model to these data.

``` r
trt1 <- rstanarm::stan_glm(log1p(mutations) ~ treatment,
                           data = md_primary_solid,
                           seed = stan_seed
                           )
trt1
```

    ## stan_glm
    ##  family:  gaussian [identity]
    ##  formula: log1p(mutations) ~ treatment
    ## ------
    ## 
    ## Estimates:
    ##                          Median MAD_SD
    ## (Intercept)              8.7    0.2   
    ## treatmenttreatment naive 0.2    0.2   
    ## sigma                    0.5    0.0   
    ## 
    ## Sample avg. posterior predictive 
    ## distribution of y (X = xbar):
    ##          Median MAD_SD
    ## mean_PPD 8.8    0.1   
    ## 
    ## ------
    ## For info on the priors used see help('prior_summary.stanreg').

This suggests the treatment effect on number of mutations may be relatively modest, with a median effect indicating that the average mutation count among treatment naive samples would be 20% higher than that among chemo-treated samples (with a relatively wide posterior interval).

``` r
bayesplot::mcmc_areas(as.array(trt1), pars = 'treatmenttreatment naive')
```

![](Hierarchical_model_mutations_and_peptides_files/figure-markdown_github-ascii_identifiers/psolid-trt-effect-1.png)

How well do this model's predictions match our data?

``` r
trt1.ppred <- rstanarm::predictive_interval(trt1) %>%
  tbl_df(.)
trt1.median <- rstanarm::predictive_interval(trt1, 0.01) %>%
  tbl_df(.) %>%
  dplyr::mutate(median = (`49.5%` + `50.5%`)/2) %>%
  dplyr::select(median)

md_primary_solid2 <- 
  md_primary_solid %>% 
  dplyr::bind_cols(trt1.ppred) %>%
  dplyr::bind_cols(trt1.median)
```

``` r
ggplot(md_primary_solid2, aes(x = treatment, y = log1p(mutations))) + 
  geom_jitter() +
  geom_errorbar(aes(x = treatment, ymin = `5%`, ymax = `95%`),
                data = md_primary_solid2 %>% dplyr::distinct(treatment, .keep_all=T),
                colour = 'red', alpha = 0.5)
```

![](Hierarchical_model_mutations_and_peptides_files/figure-markdown_github-ascii_identifiers/psolid-trt-ppred-1.png)

How well does our model recover the observed distributions of variables?

``` r
bayesplot::pp_check(trt1)
```

![](Hierarchical_model_mutations_and_peptides_files/figure-markdown_github-ascii_identifiers/psolid-glm-ppcheck-1.png)

Not bad ..

### Try a negative-binomial model?

What if we tried a negative-binomial model instead?

``` r
trt1nb <- rstanarm::stan_glm(mutations ~ treatment,
                             data = md_primary_solid,
                             family = neg_binomial_2(),
                             seed = stan_seed
)
trt1nb
```

    ## stan_glm
    ##  family:  neg_binomial_2 [log]
    ##  formula: mutations ~ treatment
    ## ------
    ## 
    ## Estimates:
    ##                          Median MAD_SD
    ## (Intercept)              8.8    0.2   
    ## treatmenttreatment naive 0.2    0.2   
    ## reciprocal_dispersion    4.3    0.7   
    ## 
    ## Sample avg. posterior predictive 
    ## distribution of y (X = xbar):
    ##          Median MAD_SD
    ## mean_PPD 7777.7  612.2
    ## 
    ## ------
    ## For info on the priors used see help('prior_summary.stanreg').

(notice that here we have almost identical parameter estimates)

``` r
bayesplot::pp_check(trt1nb)
```

![](Hierarchical_model_mutations_and_peptides_files/figure-markdown_github-ascii_identifiers/psolid-glmnb-ppcheck-1.png)

Here we have a slightly better fit, but not by much. Consistent with theory, the log-transform works well as an approximation to the 'counting process' at high levels of the counts.

### Adjust for number of cycles?

Next we look at estimating the effects of number of cycles on mutation count. Sometimes adding more information can address noise in the model, and sometimes it just .. adds noise.

``` r
trt2 <- rstanarm::stan_glm(log1p(mutations) ~ treatment + `total cycles`,
                           data = md_primary_solid %>%
                             dplyr::mutate(no_treatment = ifelse(treatment == 'treatment naive', 1, 0),
                                           treatment = ifelse(treatment != 'treatment naive', 1, 0))
                           )
trt2
```

    ## stan_glm
    ##  family:  gaussian [identity]
    ##  formula: log1p(mutations) ~ treatment + `total cycles`
    ## ------
    ## 
    ## Estimates:
    ##                Median MAD_SD
    ## (Intercept)     8.9    0.1  
    ## treatment       1.0    0.8  
    ## `total cycles` -0.2    0.1  
    ## sigma           0.5    0.0  
    ## 
    ## Sample avg. posterior predictive 
    ## distribution of y (X = xbar):
    ##          Median MAD_SD
    ## mean_PPD 8.8    0.1   
    ## 
    ## ------
    ## For info on the priors used see help('prior_summary.stanreg').

Let's plot the distributions around these effects

``` r
bayesplot::mcmc_areas(as.array(trt2), pars = c('treatment', '`total cycles`'))
```

![](Hierarchical_model_mutations_and_peptides_files/figure-markdown_github-ascii_identifiers/psolid-mvtrt-coefs-1.png)

The interpretation of these results would be that :

1.  Samples that received treatment have higher average mutation count than samples that are treatment-naive
2.  Among those receiving treatment, those with more cycles tended to have lower mutation count

This may or may not make biological sense (to me it feels like a stretch), and the posterior distributions of effects are all pretty broad. So my inclination would be to judge these effects as being "within the noise".

### Use a varying-intercept effect?

Before moving on to include solid/relapse/treated & ascites samples, we fit the same model using a varying-coefficient structure. We do this because this is the type of model we will fit in later iterations and so it will be helpful to have a baseline here to compare against.

This model is structurally very similar to fitting a no-intercept model, with formula `mutations ~ 0 + treatment`, which would estimate separate mean numbers of mutations among treated & naive solid primary samples. The difference is that, in this formulation, the treatment-specific means are drawn from a higher-level distribution of means which acts like a prior for the group-specific means, and which in effect regularizes the treatment & non-treatment means -- the group-specific means end up shrinking towards an overall inter-group mean.

``` r
trt3 <- rstanarm::stan_glmer(log1p(mutations) ~ (1 | treatment),
                           data = md_primary_solid,
                           adapt_delta = 0.999,
                           seed = stan_seed
                           )
trt3
```

    ## stan_glmer
    ##  family:  gaussian [identity]
    ##  formula: log1p(mutations) ~ (1 | treatment)
    ## ------
    ## 
    ## Estimates:
    ##             Median MAD_SD
    ## (Intercept) 8.8    0.2   
    ## sigma       0.5    0.0   
    ## 
    ## Error terms:
    ##  Groups    Name        Std.Dev.
    ##  treatment (Intercept) 0.39    
    ##  Residual              0.48    
    ## Num. levels: treatment 2 
    ## 
    ## Sample avg. posterior predictive 
    ## distribution of y (X = xbar):
    ##          Median MAD_SD
    ## mean_PPD 8.8    0.1   
    ## 
    ## ------
    ## For info on the priors used see help('prior_summary.stanreg').

The summary for this model does not include the group-specific means by default. We can, however, recover them quite easily.

``` r
bayesplot::mcmc_areas(as.array(trt3), regex_pars = '\\(Intercept\\) treatment\\:')
```

![](Hierarchical_model_mutations_and_peptides_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-1-1.png)

A textual summary can be accessed via `summary`:

``` r
summary(trt3)
```

    ## 
    ## Model Info:
    ## 
    ##  function:  stan_glmer
    ##  family:    gaussian [identity]
    ##  formula:   log1p(mutations) ~ (1 | treatment)
    ##  algorithm: sampling
    ##  priors:    see help('prior_summary')
    ##  sample:    4000 (posterior sample size)
    ##  num obs:   80
    ##  groups:    treatment (2)
    ## 
    ## Estimates:
    ##                                            mean   sd    2.5%   25%   50%
    ## (Intercept)                                8.8    0.3   8.1    8.7   8.8
    ## b[(Intercept) treatment:chemo_treated]     0.0    0.3  -0.7   -0.1   0.0
    ## b[(Intercept) treatment:treatment_naive]   0.1    0.3  -0.5    0.0   0.0
    ## sigma                                      0.5    0.0   0.4    0.5   0.5
    ## Sigma[treatment:(Intercept),(Intercept)]   0.2    0.4   0.0    0.0   0.0
    ## mean_PPD                                   8.8    0.1   8.7    8.8   8.8
    ## log-posterior                            -67.2    2.0 -72.2  -68.2 -66.8
    ##                                            75%   97.5%
    ## (Intercept)                                8.9   9.3  
    ## b[(Intercept) treatment:chemo_treated]     0.0   0.6  
    ## b[(Intercept) treatment:treatment_naive]   0.2   0.8  
    ## sigma                                      0.5   0.6  
    ## Sigma[treatment:(Intercept),(Intercept)]   0.1   1.1  
    ## mean_PPD                                   8.9   9.0  
    ## log-posterior                            -65.7 -64.4  
    ## 
    ## Diagnostics:
    ##                                          mcse Rhat n_eff
    ## (Intercept)                              0.0  1.0   766 
    ## b[(Intercept) treatment:chemo_treated]   0.0  1.0   891 
    ## b[(Intercept) treatment:treatment_naive] 0.0  1.0   766 
    ## sigma                                    0.0  1.0  2635 
    ## Sigma[treatment:(Intercept),(Intercept)] 0.0  1.0  1473 
    ## mean_PPD                                 0.0  1.0  4000 
    ## log-posterior                            0.1  1.0  1220 
    ## 
    ## For each parameter, mcse is Monte Carlo standard error, n_eff is a crude measure of effective sample size, and Rhat is the potential scale reduction factor on split chains (at convergence Rhat=1).

It's not surprising, here, that our effect has shrunken -- the mean `log1p(mutation)` count among treatment\_naive samples is 0.1 instead of 0.2. However the direction of effect is similar. Folks may reasonably disagree about which of these two results is more "correct", but the findings aren't inconsistent.

At this point, we are ready to move on to a variation of this model that includes all solid samples, estimating effects of relapse & treatment separately.

All solid tumor samples
-----------------------

Now we start to analyze a dataset including primary/untreated, primary/treated & relapse/treated samples.

First we note that the maximum number of samples per donor in this subset of our data is , meaning we have a handful of duplicate samples per donor. We will adjust for this in our analysis.

### Using a standard Bayesian glm

First we fit a standard `glm` without any donor-specific adjustments.

``` r
strt1 <- rstanarm::stan_glm(log1p(mutations) ~ treatment + timepoint,
                           data = md_solid, 
                           adapt_delta = 0.999,
                           seed = stan_seed
                           )
strt1
```

    ## stan_glm
    ##  family:  gaussian [identity]
    ##  formula: log1p(mutations) ~ treatment + timepoint
    ## ------
    ## 
    ## Estimates:
    ##                          Median MAD_SD
    ## (Intercept)              8.7    0.2   
    ## treatmenttreatment naive 0.2    0.2   
    ## timepointrecurrence      0.8    0.3   
    ## sigma                    0.5    0.0   
    ## 
    ## Sample avg. posterior predictive 
    ## distribution of y (X = xbar):
    ##          Median MAD_SD
    ## mean_PPD 8.9    0.1   
    ## 
    ## ------
    ## For info on the priors used see help('prior_summary.stanreg').

Here, we see a similar treatment effect as in our earlier analysis (which is encouraging), but the estimated "recurrence" effect is somewhat higher than the treatment effect.

How well are we recovering the distribution of our `log1p(mutation)` count in this sample?

``` r
bayesplot::pp_check(strt1)
```

![](Hierarchical_model_mutations_and_peptides_files/figure-markdown_github-ascii_identifiers/allsolid-strt1-ppcheck-1.png)

Pretty well.

Let's fit this model with an adjustment for within-id similarity.

``` r
strt2 <- rstanarm::stan_glmer(log1p(mutations) ~ treatment + timepoint + (1 | donor),
                           data = md_solid, 
                           adapt_delta = 0.999,
                           iter = 5000,
                           seed = stan_seed
                           )
```

    ## Warning: There were 4 chains where the estimated Bayesian Fraction of Missing Information was low. See
    ## http://mc-stan.org/misc/warnings.html#bfmi-low

    ## Warning: Examine the pairs() plot to diagnose sampling problems

``` r
strt2
```

    ## stan_glmer
    ##  family:  gaussian [identity]
    ##  formula: log1p(mutations) ~ treatment + timepoint + (1 | donor)
    ## ------
    ## 
    ## Estimates:
    ##                          Median MAD_SD
    ## (Intercept)              8.7    0.2   
    ## treatmenttreatment naive 0.2    0.2   
    ## timepointrecurrence      0.5    0.3   
    ## sigma                    0.3    0.1   
    ## 
    ## Error terms:
    ##  Groups   Name        Std.Dev.
    ##  donor    (Intercept) 0.36    
    ##  Residual             0.29    
    ## Num. levels: donor 81 
    ## 
    ## Sample avg. posterior predictive 
    ## distribution of y (X = xbar):
    ##          Median MAD_SD
    ## mean_PPD 8.9    0.0   
    ## 
    ## ------
    ## For info on the priors used see help('prior_summary.stanreg').

These findings are similar to those from the earlier model which did not adjust for duplicate samples per donor, although the effect is somewhat attenuated.

``` r
bayesplot::mcmc_areas(as.array(strt2), regex_pars = c('^treatment', '^timepoint'))
```

![](Hierarchical_model_mutations_and_peptides_files/figure-markdown_github-ascii_identifiers/allsolid-strt2-coefplot-1.png)

According to this model, we have an estimated average 50% increase in mutations for relapse samples, vs primary.

The 95% credible intervals for this increase are:

``` r
rstanarm::posterior_interval(strt2, prob = 0.95, pars = c('timepointrecurrence'))
```

    ##                            2.5%    97.5%
    ## timepointrecurrence -0.04411578 1.134112

How much of this probability mass is &lt; 0?

``` r
mean(sapply(as.array(strt2)[,,'timepointrecurrence'], FUN = function(x) x <= 0))
```

    ## [1] 0.0337

This number reflects the so-called 'bayesian p-value', IE the posterior probability of a relapse effect being &lt;= 0. Note that this would correspond to a one-sided p-value in traditional frequentist NHT framework.

*Side note: this adjustment for `donor` only accounts for the fact that we'd expect two samples from the same donor to be more similar from one another than two samples from different donors. It does not estimate "varying-coefficients", ie relapse or treatment effects cannot vary by donor. All effects are estimated as population averages.*

### Using a varying-coefficient model

``` r
strt3 <- rstanarm::stan_glmer(log1p(mutations) ~ treatment + timepoint + (1 + timepoint | donor),
                           data = md_solid, 
                           adapt_delta = 0.999,
                           iter = 5000,
                           seed = stan_seed
                           )
```

    ## Warning: There were 34 divergent transitions after warmup. Increasing adapt_delta above 0.999 may help. See
    ## http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup

    ## Warning: There were 4 chains where the estimated Bayesian Fraction of Missing Information was low. See
    ## http://mc-stan.org/misc/warnings.html#bfmi-low

    ## Warning: Examine the pairs() plot to diagnose sampling problems

``` r
strt3
```

    ## stan_glmer
    ##  family:  gaussian [identity]
    ##  formula: log1p(mutations) ~ treatment + timepoint + (1 + timepoint | donor)
    ## ------
    ## 
    ## Estimates:
    ##                          Median MAD_SD
    ## (Intercept)              8.7    0.2   
    ## treatmenttreatment naive 0.2    0.2   
    ## timepointrecurrence      0.5    0.4   
    ## sigma                    0.3    0.1   
    ## 
    ## Error terms:
    ##  Groups   Name                Std.Dev. Corr 
    ##  donor    (Intercept)         0.36          
    ##           timepointrecurrence 0.31     -0.08
    ##  Residual                     0.29          
    ## Num. levels: donor 81 
    ## 
    ## Sample avg. posterior predictive 
    ## distribution of y (X = xbar):
    ##          Median MAD_SD
    ## mean_PPD 8.9    0.0   
    ## 
    ## ------
    ## For info on the priors used see help('prior_summary.stanreg').

As in our previous examples, we can plot this treatment effect

``` r
bayesplot::mcmc_areas(as.array(strt3), regex_pars = c('^treatment', '^timepoint'))
```

![](Hierarchical_model_mutations_and_peptides_files/figure-markdown_github-ascii_identifiers/allsolid-plotcoefs-1.png)
