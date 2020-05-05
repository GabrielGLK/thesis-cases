# Basilisk solver installation #
git tip: *git config --global credential.helper store* -> not need input account and password everytime

<a href = "http://basilisk.fr/src/INSTALL" > Install Basilisk solver </a>

<a href = "http://basilisk.fr/src/view.h" > Install Basilisk view </a>

In addition, if you want to use Gfsview from Gerris, please use this command:
sudo apt-get install gfsview
or install from <a href = "http://gfs.sourceforge.net/wiki/index.php/Installing_from_source#Compiling_and_installing_GfsView" > source </a>:

# **Numerical simulation of single spherical bubble rising**
The <a href ="https://github.com/GabrielGLK/thesis-cases/tree/master/Basilisk-test-cases/single-bubble">single bubble </a> test cases have been proposed by <a href ="https://github.com/GabrielGLK/thesis-cases/blob/master/Basilisk-test-cases/single-bubble/reference-papers/BenchmarkBubble_IJNF2009.pdf">Hysing </a>, which are classical and typical to study and validate the accuracy of the two-phase fluid solvers. The initial configuration and physical parameters refers to read <a href = "http://www.featflow.de/en/benchmarks/cfdbenchmarking/bubble.html" > Featflow </a>.
- In Nabil's thesis, they used MooNMD as a reference solution since it is based on Lagrangian interface method.
- The <a href = "https://github.com/GabrielGLK/thesis-cases/blob/master/Basilisk-test-cases/single-bubble/2D/case-1/rising.c" > case-1 </a> with a low density and viscosity ratios, and higher surface tension coefficient which will keep the bubble more circular. 

- The <a href = "https://github.com/GabrielGLK/thesis-cases/blob/master/Basilisk-test-cases/single-bubble/2D/case-2/rising.c" > case-2 </a> with a larger density and viscosity ratios, and lower surface tension coefficient which will change the bubble shape especially at bubble tail. 


# **Numerical simulation of Taylor bubble rising in stagnant liquid** 
In this part, we introduce the test cases using Basilisk to simulate Taylor bubble rising in quiescent liquid through singularities. The comparisons have been done between experiments results from Chengsi and simulations results, including bubble rising in straight tube and expansion/contraction tube.

In addtion, we extend the simulation cases into several new singularities, including gradual, parabolic concave/convex expansions/contractions. I will introduce more in the following parts:

