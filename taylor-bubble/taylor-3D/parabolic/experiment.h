// experimental data comaring with chengsi

#define percent_60 0
#define percent_70 0
#define percent_80 1
#define percent_90 0
#define percent_95 0


#define case_1 0
#define case_2 1
#define case_3 0
#define case_4 0
#define case_5 0
#define case_6 0
#define case_7 0
#define case_8 0

// different change-ratio
#if case_1
#define ratio 1.72
#endif

#if case_2
#define ratio 1.45
#endif

#if case_3
#define ratio 1.33
#endif

#if case_4
#define ratio 1.24
#endif

#if case_5
#define ratio 1.12
#endif

#if case_6
#define ratio 0.9
#endif

#if case_7
#define ratio 0.81
#endif

#if case_8
#define ratio 0.69
#endif

/********************************************/
// fluid peroperties
/*
length_ratio(bubble tail to the bottom(4r0)+bubble body length(not including head)
This is compared with chengsi's experiment
*/


/************************* 60% *******************************/
#if percent_60
#define RHOR 910
#define MUR 1655
#define N_f 230.9
#define Eo 40.55

#if case_1
#define length_ratio 2.77 
#define r_0 0.421
#endif 
#if case_2
#define length_ratio 3.71
#define r_0 0.422
#endif 
#if case_3
#define length_ratio 4.54 
#define r_0 0.416
#endif 
#if case_4
#define length_ratio 3.54
#define r_0 0.424
#endif 
#if case_5
#define length_ratio 3.88
#define r_0 0.422
#endif 
#if case_6
#define length_ratio 2.68 
#define r_0 0.42
#endif 
#if case_7
#define length_ratio 3.77
#define r_0 0.424
#endif 
#if case_8
#define length_ratio 2.66
#define r_0 0.417
#endif 

#endif

/************************* 70% *******************************/

#if percent_70
#define RHOR 934
#define MUR 2906
#define N_f 134.93
#define Eo 41.86

#if case_1
#define length_ratio 3.25
#define r_0 0.4
#endif 
#if case_2
#define length_ratio 3.47
#define r_0 0.4
#endif 
#if case_3
#define length_ratio 4.73
#define r_0 0.4
#endif 
#if case_4
#define length_ratio 2.94
#define r_0 0.404
#endif 
#if case_5
#define length_ratio 3.99 
#define r_0 0.4
#endif 
#if case_6
#define length_ratio 6.02
#define r_0 0.395
#endif 
#if case_7
#define length_ratio 5.23 
#define r_0 0.405
#endif 
#if case_8
#define length_ratio 2.21 
#define r_0 0.403
#endif 

#endif

/************************* 80% *******************************/

#if percent_80
#define RHOR 1044
#define MUR 5350
#define N_f 81.91
#define Eo 47.52

#if case_1
#define length_ratio 3.9
#define r_0 0.38
#endif 
#if case_2
#define length_ratio 4.04
#define r_0 0.383
#endif 
#if case_3
#define length_ratio 4.78
#define r_0 0.37
#endif 
#if case_4
#define length_ratio 6.37
#define r_0 0.386
#endif 
#if case_5
#define length_ratio 4.32 
#define r_0 0.38
#endif 
#if case_6
#define length_ratio 3.6
#define r_0 0.38
#endif 
#if case_7
#define length_ratio 3.75 
#define r_0 0.378
#endif 
#if case_8
#define length_ratio 3.9 
#define r_0 0.379
#endif 
#endif

/************************* 90% *******************************/

#if percent_90
#define RHOR 1067
#define MUR 16344
#define N_f 27.4
#define Eo 49.1
#define Mo 0.21

#if case_1
#define length_ratio 7.4
#define r_0 0.354
#endif 
#if case_2
#define length_ratio 5.3
#define r_0 0.354
#endif 
#if case_3
#define length_ratio 5.2
#define r_0 0.347
#endif 
#if case_4
#define length_ratio 6
#define r_0 0.352
#endif 
#if case_5
#define length_ratio 8
#define r_0 0.352
#endif 
#if case_6
#define length_ratio 4
#define r_0 0.351
#endif 
#if case_7
#define length_ratio 3.8
#define r_0 0.348 
#endif 
#if case_8
#define length_ratio 3.6
#define r_0 0.345 
#endif 
#endif

/************************* 95% *******************************/

#if percent_95
#define RHOR 1075
#define MUR 29809
#define N_f 15.1
#define Eo 49.8

#if case_1
#define length_ratio 4.5
#define r_0 0.341
#endif 
#if case_2
#define length_ratio 5.67
#define r_0 0.341
#endif 
#if case_3
#define length_ratio 5
#define r_0 0.347
#endif 
#if case_4
#define length_ratio 3.8
#define r_0 0.364
#endif 
#if case_5
#define length_ratio 5.84
#define r_0 0.344 
#endif 
#if case_6
#define length_ratio 5
#define r_0 0.343
#endif 
#if case_7
#define length_ratio 6.83
#define r_0 0.352
#endif 
#if case_8
#define length_ratio 3
#define r_0 0.361

#endif 
#endif 
