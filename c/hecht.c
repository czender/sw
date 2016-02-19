/* Charlie,
	I found a simple program I used to compute loan schedules.  I
set it up to compute schedules for either fixed rate or arm (using the
worst case scenario for the arm, going up the maximum amount each
year).  I had a term of 180 months (15 years) which you may want to
change to 360 months using your editor, and since most 30 year arms
can reset by 2 points per year you might want to change "rate += 0.01"
to "rate += 0.02".
	Here it is:

-----Matthew.
-- 
Matthew Hecht
National Center for Atmospheric Research
PO Box 3000				phone: (303) 497-1702
Boulder, CO 80307-3000                  e-mail: hecht@ncar.ucar.edu
*/

#include <stdio.h>
#include <math.h>
int main()
{
  float principal, rate, rate_max, payment, twelfth, pay();
  float nbr_year;
  float interest;
  float tot_int=0.0;
  int n_months, count, arm_flag, year;

  twelfth=1.0/12.0;

  printf("enter principal:");
  scanf("%f", &principal);
  printf("enter initial interest rate (0.07, etc):");
  scanf("%f", &rate);
  printf("enter # years (15, etc):");
  scanf("%f", &nbr_year);
  printf("enter 0 for conventional, 1 for 1/6 arm:");
  scanf("%d", &arm_flag);

  rate_max = rate + 0.06;
  n_months = (int)nbr_year*12;
  year = 0;
  payment = pay(principal, rate, n_months);

  while (n_months > 0) {
    if (arm_flag == 1 && year > 0 && rate < 0.99*rate_max){
      rate += 0.01;
      payment = pay(principal, rate, n_months);
    }
    
    for (count=1; count<=12; count++, n_months--){
      interest = rate*principal*twelfth;
      tot_int += interest;
      principal -= (payment-interest);
    }
    year++;
    printf("after year %2d, at %f,\n",
	   year, rate);
    printf("\t payment = %6.2f, principal = %6.2f, interest = %6.2f, total principal = %6.2f, total interest = %6.2f\n",
	   payment, payment-interest, interest, principal, tot_int);
  }
  return 0;
}

float pay(principal, rate, n_months)
     float principal, rate;
     int n_months;
{
  float payment, twelfth, interest, tot_int, residual;
  float abs_res, f_accel;
  float resid();
  int n_m, arm_flag, it_count;
  
  twelfth=1.0/12.0;
  f_accel = 0.5;
  n_m = 0;
  tot_int = 0.0;
  
  /* a starting point... */
  payment = 1.25*rate*principal*twelfth;
  /* printf("using payment of %f\n", payment); */
  
  residual = resid( n_months, principal, rate, payment );
  if (residual > 0.0) {
    abs_res = residual;
  } else {
    abs_res = -residual;
  }

  it_count = 1;
  while (abs_res > 10.0) {
    it_count++;
    payment += f_accel * residual / (float)n_months;
    /* printf("new payment is %f...\n", payment); */
    
    residual = resid( n_months, principal, rate, payment );
    if (residual > 0.0) {
      abs_res = residual;
    } else {
      abs_res = -residual;
    }
  }
  /* printf("after %d iterations, remaining principal is %f.\n",
     it_count, residual);*/
  return payment;
}

float resid( n_months, principal, rate, payment )
     int n_months;
     float principal, rate, payment;
{
  int n_m;
  float twelfth, interest;
  twelfth=1.0/12.0;

  for (n_m=1; n_m<=n_months; n_m++) {
    if (principal > 0.0) {
      interest = rate*principal*twelfth;
    } else {
      interest = 0.0;
    }
    principal -= (payment-interest);
  }

  return principal;
}
