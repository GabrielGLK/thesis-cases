// experimental data comaring with chengsi

#define percent_60 1
#define percent_70 0
#define percent_80 0
#define percent_90 0
#define percent_95 0


#define case_1 1
#define case_2 0
#define case_3 0
#define case_4 0
#define case_5 0
#define case_6 0
#define case_7 0
#define case_8 0
#define case_9 0
// different change-ratio
#if case_1
//32mm/28mm 20mm/16.2mm
#define ratio 1.73
#endif

#if case_2
//27.9mm/23.9mm 20mm/16.2mm
#define ratio 1.475
#endif

#if case_3
//26mm/21.95mm 20mm/16.2mm
#define ratio 1.355
#endif

#if case_4
//23.85mm/20.15mm 20mm/16.2mm
#define ratio 1.244
#endif

#if case_5
//21.87mm/18.2mm 20mm/16.2mm
#define ratio 1.123
#endif

#if case_6
//18.9mm/15.2mm 20mm/16.2mm
#define ratio 0.94
#endif

#if case_7
//17.15mm/13.45mm 20mm/16.2mm
#define ratio 0.83
#endif

#if case_8
//15.05mm/11.30mm 20mm/16.2mm
#define ratio 0.7
#endif

// straght tubes
#if case_9
#define ratio 1
#endif
/********************************************/
// fluid peroperties
/*
length_ratio(bubble tail to the bottom(4r0)+bubble body length(not including head)
This is compared with chengsi's experiment
*/


/************************* 60% *******************************/
#if percent_60
#define RHOR 997.0 
#define MUR 1651.4
#define N_f 261.3 
#define Eo 45.1
#define Mo 1.96e-5

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

#if case_9
#define length_ratio 3.77
#define r_0 0.424
#endif 

#endif

/************************* 70% *******************************/

#if percent_70
#define RHOR 1018.4 
#define MUR 2900.6
#define N_f 152.0
#define Eo 46.7
#define Mo 0.000191

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

#if case_9
#define length_ratio 3.99 
#define r_0 0.4
#endif 

#endif

/************************* 80% *******************************/

#if percent_80
#define RHOR 1038.8
#define MUR 5340.6
#define N_f 84.2
#define Eo 48.3
#define Mo 0.00224

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
#if case_9
#define length_ratio 3.75 
#define r_0 0.378
#endif 

#endif

/************************* 90% *******************************/

#if percent_90
#define RHOR 1058.2
#define MUR 16315.4
#define N_f 27.4
#define Eo 49.8
#define Mo 0.199

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
#if case_9
#define length_ratio 3.8
#define r_0 0.348 
#endif
#endif

/************************* 95% *******************************/

#if percent_95
#define RHOR 1067.4
#define MUR 29757.7
#define N_f 15.5
#define Eo 50.6
#define Mo 2.23

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
#if case_9
#define length_ratio 3.8
#define r_0 0.364
#endif
#endif 
