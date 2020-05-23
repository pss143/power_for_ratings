library(MASS)
library(lme4)
library(lmerTest)

#' Simulate a study in which subjects rate items for a trait, and
#'  the items vary in their actual value for this trait. The true correlation
#'  between item actual values and subject ratings is systematically
#'  varied, and we assess the power of different designs to detect
#'  correlations that are present
#'  Based on DeBruine & Barr (2019)
#'  https://debruine.github.io/lmem_sim/index.html
#'
#' NB A subject never rates an item more than once.
#' This simplifies everything

#' Subject ratings are generated and assessed by the lme model below.
#' These models are ideal when we have properly two random variables.
#'   rating = b0 + I0i + S0s + ((b1 + S1s) * actual) + err, where:
#'.    rating = rating for item i made by subject s
#'.    b0 = grand mean rating (== 0)
#'.    I0i = random intercept for item i
#'     S0s = random intercept for subject s
#'     b1 = the true correlation of item value and subject ratings
#'     S1s = the random slope for subject s, their sensivity to any
#'       true signal in items
#'.    actual = the actual value for item i
#'.    err = residual error
#' These paramters are fixed a priori to produce fixed and random effects
#'  with normal distributions of mean=0, sd=1
#' The exception is b1, the true correlation of actual and ratings,
#'  which is manipulated in the sims to test our ability to detect small,
#'  medium, large effect sizes

###### Manipulated variables
#' Given the data generating process as fixed, the variables below are what we
#' can manipulate to optimise the power of the study. These are passed to the
#' function full_sim()
#'
#' The hard upper limit on the power of our study is the number of stimulus
#' items. This can't readily be increased, and so we intend to use our full set
#' of 1000 items, and so this value is fixed in the simulations.
nitems = 1000

#' We assess power on "true" correlations of different sizes, including 0. The
#' correlation of actual item value and subject ratings, ie real and perceived
#' value. How frequent are false alarms when r=0, how low can r go and still be
#' reliably detected?
true_r = c(0, .5)

#' Given all this, there are only 2 design variables that can affect
#' experimental power:
#' trials_per_subject and ratings_per_item
#'
#' trials_per_subject (SCALAR). Since a subject never rates an item twice, this
#' can also be seen as the number of items rated per subject. For example, if
#' trials_per_subject = 1, each subject rates just a single item; if
#' trials_per_subject = nitems, then the design is completely within-subject, as
#' each subject will rate all items trials_per_subject can be set to reflect
#' what we can realistically ask subjects to do, eg in an online vs lab-based
#' study. For online, we want this to be dozens not 100s.
trials_per_subject = 40

#' ratings_per_item (VECTOR). The mean rating for each item will consist of this
#' many values, each from a different subject, with their own random intercepts
#' and slopes. Obviously, as the number of ratings increases, the mean rating of
#' an item will approach its actual value. The big question is, how many ratings
#' are needed. Specified as a vector to test a range of values
ratings_per_item = c(5, 25, 50)

#' Simulation control
sims_per_cond = 1 #' The number of sims per condition.
#' The number of conditions equals length(ratings_per_item) * length(true_r)
verbose = TRUE #' diagnostics or no

#' run one simulated experiment
single_expt = function(nitems, ratings_per_item, trials_per_subject,
                       true_r, verbose = FALSE) {
  #' data-generation parmeters -- fixed
  actual_mu = 0 # mean of actual item values (subject trying to estimate these)
  actual_sd = 1 # sd of item values
  b0 = 0 # mean rating, boring
  I0i_sd = 1 # sd of by-item random intercepts (mean = 0)
  S0s_sd = 1 # sd of by-subject random intercepts (mean = 0)
  S1s_sd = 1 # sd of by-subject random slopes of ratings with actual
  S01_cor = 0
  #' S01_cor, the correlation of by-subject intercepts and slopes is set
  #' to 0. Although this correlation can be interesting in some contexts, it
  #'  doesn't make much sense here. No prior reason to
  #'  expect observers who are more sensitive to actual item value to have
  #'  a systematic high or low bias in their overall ratings. But **just in case** someone wants to model this,
  #'  the code in sample_subjects uses the approach from DeBruine and Barr
  err_sd = 1 # sd of residual error

  #' sampleitems -- Simulate the sampling processes for items
  #' Randomly assign actual values for each item, as well as
  #'  the item-specific random intercept I0i
  sample_items = function(nitems) {
    items = data.frame(id = 1:nitems,
      actual = rnorm(n = nitems, mean = actual_mu, sd = actual_sd),
      I0i = rnorm(n = nitems, mean = 0, sd = I0i_sd)
    )
    return(items)
  }

  #' samplesubjects -- Simulate sampling of subjects, specifically, generate
  #'  the subject-specific intercepts and slopes
  sample_subjects = function(nsubjects) {
    #' The only complication below is allowing for the by-subject random intercepts
    #' andvslopes to be correlated. This might be generally good practice, but for
    #' this specific case, see the notes above on the S01_cor parameter. But **just in
    #' case** someone wants to model this, the code below allows it using the
    #' corr_mat and var_mat approach taken from DeBruine & Barr
    corr_mat = matrix(c(1, S01_cor, S01_cor, 1), nrow = 2, byrow = TRUE)
    # (multiply the SDs for each cell, so SD of each column reflect S0/1_sd)
    var_mat <- matrix(c(S0s_sd * S0s_sd, S0s_sd * S1s_sd,
                        S0s_sd * S1s_sd, S1s_sd * S1s_sd),
                      nrow = 2, byrow = TRUE)
    S <- MASS::mvrnorm(n = nsubjects, mu = c(0, 0),
                       Sigma = corr_mat * var_mat, empirical = FALSE)
    subjects = data.frame(id = 1:nsubjects, S0s = S[, 1], S1s = S[, 2])
    return(subjects)
  }

  ntrials = nitems * ratings_per_item
  nsubjects = ntrials / trials_per_subject

  if(verbose) {
    outstr = paste('true r = %.2f',
                   'nitems: %d; trials per item): %d',
                   'nsubjects: %d; trials per subject: %d',
                   'total trials: %d', sep = '\n')
    message(sprintf(outstr, true_r, nitems, ratings_per_item,
                    nsubjects, trials_per_subject, ntrials))
  }
  items = sample_items(nitems)
  subjects = sample_subjects(nsubjects)
  trials = data.frame(
    item = rep(1:nitems, times = ratings_per_item),
    subject = rep(1:nsubjects, each = trials_per_subject),
    b0 = b0,
    b1 = true_r,
    err = rnorm(n = ntrials, mean = 0, sd = err_sd))
  # merge with the random effects generated in sample_items and sample_subjects
  trials = merge(trials, items, by.x = 'item', by.y = 'id')
  trials = merge(trials, subjects, by.x = 'subject', by.y = 'id')
  trials = trials[, c('subject', 'item', 'b0' ,'I0i', 'S0s', 'b1', 'S1s', 'actual', 'err')]

  trials$rating = with(trials,
                       b0 + I0i + S0s + ((b1 + S1s) * actual) + err)
  return(trials)
}

full_sim = function(nitems, ratings_per_item, trials_per_subject,
                    true_r, sims_per_cond, verbose = FALSE) {
  if((nitems*ratings_per_item) %% trials_per_subject != 0) {
    stop('Design uses fractional number of subjects')
  }
  #' pre-define simulation data structure to hopefully speed things along
  nsims = length(ratings_per_item) * length(true_r) * sims_per_cond
  sim_results = data.frame(nitems = rep(nitems, nsims),
                           nsubjects = NA,
                           ratings_per_item = NA, trials_per_subject = NA,
                           true_r = NA, beta = NA, se = NA,
                           t = NA, df = NA, p = NA)
  start_time = Sys.time()
  sims_so_far = 0
  for(b1 in true_r) {
    for(rpi in ratings_per_item) {
      ntrials = nitems * rpi
      nsubjects = ntrials / trials_per_subject
      for(sim in 1:sims_per_cond) {
        sims_so_far = sims_so_far + 1
        if(verbose){ message(sprintf('sim %d of %d', sims_so_far, nsims))}
        trials = single_expt(nitems, rpi, trials_per_subject, b1, verbose)
        #' NB As we generating vast amounts of simulated data of known
        #' characteristics, the tolerance for convergence is more liberal
        fit = lmer(rating ~ 1 + (1|item) + (actual||subject) + actual,
                  data = trials,
                  control = lmerControl(check.conv.grad =
                          .makeCC("warning", tol = 5e-2, relTol = NULL)),
                  REML=FALSE)
        out = summary(fit)$coefficients[2, ] # all the info for the fixed effect b1
        sim_results[sims_so_far, ] = c(nitems, nsubjects, rpi, trials_per_subject,
                                       b1, out)
        print(sim_results[sims_so_far, ])
      }
    }
  }
  stop_time = Sys.time()
  message('Done')
  print(stop_time - start_time)
  return(sim_results)
}

#' indicative plot of estimated betas
betaplot = function(simout, alpha=.01) {
  true_rs = as.numeric(attr(table(simout$true_r), 'dimnames')[[1]])
  legtext = as.character(true_rs)
  rdat = subset(simout, true_r == true_rs[1])
  plot(density(rdat$beta),
       xlim = c(min(simout$beta), max(simout$beta)),
       main = 'Estimated beta as a function of true r')
  legtext[1] = sprintf('r=%.2f detect: %.2f', true_rs[1],
                       nrow(subset(rdat, p < alpha))/nrow(rdat))
  for(i in 2:length(true_rs)) {
    rdat = subset(simout, true_r == true_rs[i])
    lines(density(rdat$beta), col=i)
    legtext[i] = sprintf('r=%.2f detect: %.2f', true_rs[i],
                         nrow(subset(rdat, p < alpha))/nrow(rdat))
  }
  legend(legend = legtext, x='topright')
  #'The criterion line drawn here is just a rough indicator of where the cutoff
  #'for significance at p=.01 would lie. It is assuming that in these sims the
  #'ns are almost inf, and se for beta estimates doesn't change much.
  abline(v = mean(simout$se)*2.56)
}

out1 = full_sim(nitems=1000, ratings_per_item = 30,
               trials_per_subject = 40,
               true_r = c(0, .1, .2),
               sims_per_cond = 100, verbose = TRUE)

#' In this plot, the criterion line is approximately at a beta estimate
#' for alpha = .01. When r=0 these are false positives
#' For other values of r the ability to detect the effect that is
#' truly present is indicated.
#' In the legend, the detection rates for each value of true_r are given
betaplot(out1)
