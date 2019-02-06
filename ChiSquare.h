#pragma once

#include <iostream>
#include <cmath>
 
double lowerincompletegamma(double shape, double x)
{
       double sum = 0;
       double term = 1.0/shape;
       double n = 1;
 
       while (term != 0 && n < 100) {
              sum=sum+term;
              term = term*(x/(shape+n));
              n++;
       }
 
       return pow(x,shape)*exp(-1*x)*sum;
}
 
double gammln(double xx) // From numerical recipes in C
{
       double x, y, tmp, ser;
       static double cof[6]={ 76.18009172947146, -86.50532032941677, 24.01409824083091,
              -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5 };
       int j;
 
       y = x = xx;
       tmp = x+5.5;
       tmp -= (x+0.5)*log(tmp);
       ser = 1.000000000190015;
       for (j = 0; j <= 5; j++) ser += cof[j]/++y;
       return -tmp+log(2.5066282746310005*ser/x);
}
 
double chisqr(double teststat, double df)
{
       double numerator;
       double denominator;
 
       numerator = lowerincompletegamma(df/2,teststat/2);
       denominator = exp(gammln(df/2));
 
       return numerator/denominator;
}
 
 