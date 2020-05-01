// difference fluid properties for comparising with reference

#define A 4.7e-5
#define B 1.6e-2
#define C 0.33
#define D 2.8
#define E 105
#define F 2.9e6

#define A_A 1
#define B_B 0
#define C_C 0
#define D_D 0
#define E_E 0
#define F_F 0

#define case1 0
#define case2 0
#define case3 0
#define case4 0
#define case5 0
#define case6 0
#define case7 0
#define case8 1
#define case9 0
#define case10 0



/************************* M0 = 4.7e-5 ********************************/
#if A_A
#if case1
#define N_f 43.4065163806
#define Eo 5.5051927529  
#endif

#if case2
#define N_f 56.9510918849
#define Eo 7.9074300579  
#endif

#if case3
#define N_f 62.971090219
#define Eo 9.0410922883 
#endif

#if case4
#define N_f 82.8684919407
#define Eo 13.0382268608 
#endif

#if case5
#define N_f 95.6884909041
#define Eo 15.7947307341
#endif

#if case6
#define N_f 173.5429783137
#define Eo 34.9334931354 
#endif

#if case7
#define N_f 401.6671606011
#define Eo 106.9517071151 
#endif

#if case8
#define N_f 111
#define Eo 187.035
#endif
#endif


/************************* M0 = 1.6e-2 ********************************/
#if B_B
#if case1
#define N_f 11.2338661792
#define Eo 6.3398500937 
#endif

#if case2
#define N_f 14.6814385399
#define Eo 9.0586777272 
#endif

#if case3
#define N_f 16.4395210034
#define Eo 10.5331646613 
#endif

#if case4
#define N_f 21.6965281274
#define Eo 15.2484951461
#endif

#if case5
#define N_f 25.2452524181
#define Eo 18.6614886111
#endif

#if case6
#define N_f 46.7135569993
#define Eo 42.3933019766
#endif

#if case7
#define N_f 71.9747449285
#define Eo 75.4420544776
#endif

#if case8
#define N_f 109.3957497045
#define Eo 131.8379959527 
#endif
#endif

/************************* M0 = 0.33 ********************************/
#if C_C
#if case1
#define N_f 7.0144616934
#define Eo 9.2789131637 
#endif

#if case2
#define N_f 7.7601774658
#define Eo 10.6169600743
#endif

#if case3
#define N_f 10.142385719
#define Eo 15.1713703648
#endif

#if case4
#define N_f 11.8233413043
#define Eo 18.6133673541
#endif

#if case5
#define N_f 21.5451528861
#define Eo 41.428965311
#endif

#if case6
#define N_f 32.7913489815
#define Eo 72.529890903
#endif

#if case7
#define N_f 49.9539427416
#define Eo 127.1347701586
#endif
#endif

/************************* M0 = 2.8 ********************************/
#if D_D
#if case1
#define N_f 3.382486977
#define Eo 7.1565029511
#endif

#if case2
#define N_f 4.2835144447
#define Eo 9.8051217678
#endif

#if case3
#define N_f 5.5323716462
#define Eo 13.7911708786
#endif

#if case4
#define N_f 6.1952187545
#define Eo 16.0371874375
#endif

#if case5
#define N_f 8.2818054102
#define Eo 23.6167087213
#endif

#if case6
#define N_f 9.5519234659
#define Eo 28.5654096027
#endif

#if case7
#define N_f 17.7550378832
#define Eo 65.2852114113
#endif

#if case8
#define N_f 26.8416236497
#define Eo 113.2742131657
#endif

#if case9
#define N_f 41.1818847303
#define Eo 200.4445728578
#endif

#if case10
#define N_f 72.8734728192
#define Eo 429.0210806524
#endif
#endif

/************************* M0 = 105 ********************************/
#if E_E
#if case1
#define N_f 2.7304041972
#define Eo 18.0038626127
#endif

#if case2
#define N_f 4.7503051485
#define Eo 37.6725580822
#endif

#if case3
#define N_f 7.4664169374
#define Eo 68.8461116832
#endif

#if case4
#define N_f 11.3760745222
#define Eo 120.7034464896
#endif
#endif

/************************* M0 = 2.9e6 ********************************/
#if F_F
#if case1
#define N_f 0.5964989509
#define Eo 71.6053675365
#endif

#if case2
#define N_f 0.9111754945
#define Eo 125.9704577466
#endif

#if case3
#define N_f 1.6210862151
#define Eo 271.5652210076
#endif
#endif