// experimental data comaring with chengsi

#define percent_60 0
#define percent_70 0
#define percent_80 1
#define percent_90 0
#define percent_95 0
#define reference 0


#define case_1 1
#define case_2 0
#define case_3 0
#define case_4 0
#define reference 0

// different change-ratio
#if case_1
#define ratio 0.95
#endif

#if case_2
#define ratio 0.9
#endif

#if case_3
#define ratio 0.81
#endif

#if case_4
#define ratio 0.69
#endif



#if reference
#define ratio 1
#endif

/********************************************/
// fluid peroperties
/*
length_ratio(bubble tail to the bottom(4r0)+bubble body length(not including head)
This is compared with chengsi's experiment
*/
double l_taylor = 6;
#define length_ratio l_taylor

/************************* 60% *******************************/
#if percent_60
double RHOR = 910;
double MUR = 1655;
double N_f =  230.9;
double Eo = 40.55;
#endif

/************************* 70% *******************************/

#if percent_70
double RHOR = 934;
double MUR = 2906;
double N_f = 134.93;
double Eo = 41.86;
#endif

/************************* 80% *******************************/

#if percent_80
double RHOR = 1044;
double MUR = 5350;
double N_f = 81.91;
double Eo = 47.52;
#endif

/************************* 90% *******************************/

#if percent_90
double RHOR = 1067;
double MUR = 16344;
double N_f = 27.4;
double Eo = 49.1;
#endif

/************************* 95% *******************************/

#if percent_95
double RHOR = 1075;
double MUR = 29809;
double N_f = 15.1;
double Eo = 49.8;
#endif 


#if reference

#define length_ratio 6
#define RHOR 1000
#define MUR 100
#define N_f 50
#define Eo 200
#endif 
