#!/bin/bash
for replicate in {1..30}
do
  R CMD BATCH --no-save --no-restore '--args '"${replicate}" MCMC.r MCMC-final-${replicate}.Rout &
done
