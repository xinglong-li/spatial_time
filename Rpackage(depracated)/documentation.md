## API

### Input:

This model fits 2 mixed effects model, one for observation and one for palcement.
Each model is a mixed effects model, so the input should be :
 - one dataframe for fixed effects of regression
 - one dataframe for random effects of regression
 
 - one dataframe for fixed effects of palcement 
 - one dataframe for random effects of palcement

## Steps in model fitting:

### 1. Formula decode/encode

Since in this package we only focus on spatio-temporal data, the formula should be simpler than that of INLA? One problem with INLA is that it only accepts one formula, even there are more than one likelihoods, which makes formula of INLA hard to read/write.

Ideally we ask users to input independent formulas for observation and site selection, which is more intuitive.

We can also make use of the fact that our model is restricted to spatio-temporal data, so maybe we can restrict some terms in the formula, e.g. "Time", "Space"...

Maybe we can mimic the formula of the lmer function in fitting mixed effects models.

Internally we should decode the formula inputted by the user and encode it to fit INLA.

To begin with we can use INLA formula format, and to see how to simplify the form to make it.

### 2. Data standardization/rescale

Standardize the data before fitting model and rescale the result after model fitting.

### 3. Build model

Build components required by INLA, e.g. mesh, projection matrix, etc.

### 4. Model fitting

Fit the model using INLA.

## Post model fitting

### - Prediction

### - Summary/Print

### - Plot

## ISSUES:

-   How to set `boundary, max.edge, cutoff, offset` in calling `inla.mesh.2d`. and

    how to set `convex` in calling `inla.nonconvex.hull.`

-   In defining the formula for the joint likelihood, the key difference from the separate model is that in the joint model the effects can be shared among likelihoods.
