# SIV
Synthetic Instrumental Variable Method
This paper introduces a novel Synthetic Instrumental Variable (SIV) method that constructs valid instruments using only existing data, addressing key challenges in traditional instrumental variable approaches. We demonstrate that any valid instrument can be represented as a linear combination of coplanar vectors spanned by the outcome and endogenous variables in a reduced form. Based on these coplanar vectors, we develop a "dual tendency" (DT) condition that provides a moment-based criterion for identifying valid SIVs, bridging the gap between unobservable orthogonality conditions and observable data characteristics. The SIV method is robust to address potential heteroscedasticity. Our method can also determine the true sign of $cov(\bf x,u)$, often assumed a priori in empirical work.  Using simulated data and empirical applications, we demonstrate the effectiveness of our approach. The SIV method offers advantages over traditional IV approaches by mitigating issues of weak or invalid instruments and reducing reliance on scarce external instruments. This approach has broad implications for improving causal inference in various fields, such as economics, epidemiology, and policy evaluation.


There two versions of the code: 1. SIV_template code.R computes only SIV estimations. 2. SIV template code with external IV.R computes also the model using your own external IVs and puts the reslts togetehr with SIV estimations for comparison.
To use the code for your own regression model, you need to update this part of the code:
####### DATA 
Here you need to upload your own dataset instead of the code given below.
library(wooldridge)
mydata<-wooldridge::mroz

############Setting the regression variables. You need to input your variables here
H0 <- data.frame(hours, lwage,educ, age,kidslt6, kidsge6, nwifeinc)  <- instead of these variables you lsit your own variables.
Vector of variables used in the regression. First var is outcome, the second var is endogenous regressor, the rest is the list of exogenous vars.
