/*  This is a recreation in C++ of Jay Forrester's World2 program.

    The original World2 was written in the DYNAMO language and appears in
    Appendix B, Equations of the World Model, in the book World Dynamics,
    by Jay W. Forrester, 1971.

    World2 models the demographic, industrial, agricultural and natural
    resource sources and sinks of the earth.

    I hereby place this code in the public domain or, if you prefer, I
    release it under either CC0 1.0 Universal or the MIT License. All
    quoted text remains the copyright of the respective copyright owner.
    Anthony Hay, 2021, Devon, UK
*/


#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <cstring>
#include <cfloat>



namespace micro_test_library {

unsigned g_test_count;      // total number of tests executed
unsigned g_fault_count;     // total number of tests that fail


// write a message to std::cout if !(value == expected_value)
// (this is a function so that value and expected_value are evaluated only once)
template<typename A, typename B>
void test_equal(const A & actual_value, const B & expected_value,
    const char * filename, const size_t linenum, const char * function_name)
{
    ++g_test_count;
    if (!(actual_value == expected_value)) {
        ++g_fault_count;
        // e.g. love.cpp(2021) : in proposal() expected 'Yes!', but got 'Hahaha'
        std::cout
            << filename
            << '('              << linenum
            << ") : in "        << function_name
            << "() expected '"  << expected_value
            << "', but got '"   << actual_value
            << "'\n";
    }
}

#define TEST_EQUAL(value, expected_value)                   \
{                                                           \
    micro_test_library::test_equal(value, expected_value,   \
        __FILE__, __LINE__, __FUNCTION__);                  \
}


void test_equal_double(const double & actual_value, const double & expected_value,
    const char * filename, const size_t linenum, const char * function_name)
{
    ++g_test_count;
    // N.B. this is *not* a general solution, but sufficient for our needs here
    if (!(std::fabs(actual_value - expected_value) < DBL_EPSILON)) {
        ++g_fault_count;
        std::cout
            << filename
            << '('              << linenum
            << ") : in "        << function_name
            << "() expected '"  << expected_value
            << "', but got '"   << actual_value
            << "'\n";
    }
}

#define TEST_EQUAL_DOUBLE(value, expected_value)                \
{                                                               \
    micro_test_library::test_equal_double(value, expected_value,\
        __FILE__, __LINE__, __FUNCTION__);                      \
}



} //namespace mico_test_library







////////  //    // //    //    ///    //     //  ///////  
//     //  //  //  ///   //   // //   ///   /// //     // 
//     //   ////   ////  //  //   //  //// //// //     // 
//     //    //    // // // //     // // /// // //     // 
//     //    //    //  //// ///////// //     // //     // 
//     //    //    //   /// //     // //     // //     // 
////////     //    //    // //     // //     //  /////// 
namespace dynamo {



// simulate DYNAMO CLIP() function
double clip(double a, double b, double c, double d)
{
    return c >= d ? a : b;
}


// simulate DYNAMO TABHL() function
// return y value for given 'x' using linear interpolation of given data set
// (xstart, ytbl[0]), (xstart+xstep, ytbl[1]) ... (xend, ytbl.back());
// use the extreme value in 'ytbl' when range exceeded; modeled on DYNAMO TABHL function
double tabhl(const std::vector<double> & ytbl, double x, double xstart, double xend, double xstep)
{
    const size_t range = static_cast<size_t>((xend - xstart) / xstep + 1);
    if (ytbl.size() != range)
        throw std::runtime_error("tabhl() size of given 'ytbl' does not match given xrange");

    if (xstart < xend) {
        if (x < xstart)
            return ytbl[0];
        if (x > xend)
            return ytbl.back();
    }
    else {
        if (x < xend)
            return ytbl.back();
        if (x > xstart)
            return ytbl[0];
    }

    const size_t i = static_cast<size_t>((x - xstart) / xstep);
    if (i == range - 1)
        return ytbl[i];

    // linear interpolation: y = y0 + (x - x0) * (y1 - y0) / (x1 - x0)
    return ytbl[i] + (x - xstart - (xstep * i)) * (ytbl[i + 1] - ytbl[i]) / xstep;
}


// simulate DYNAMO TABLE() function
// return y value for given 'x' using linear interpolation of given data set
// (xstart, ytbl[0]), (xstart+xstep, ytbl[1]) ... (xend, ytbl.back());
// requires that x lies between xstart and xend; modeled on DYNAMO TABLE function
double table(const std::vector<double> & ytbl, double x, double xstart, double xend, double xstep)
{
    if (xstart < xend) {
        if (x < xstart || x > xend)
            throw std::runtime_error("table() given 'x' out of range");
    }
    else {
        if (x < xend || x > xstart)
            throw std::runtime_error("table() given 'x' out of range");
    }
    return tabhl(ytbl, x, xstart, xend, xstep);
}

}//namespace dynamo



//      //  ///////  ////////  //       ////////   ///////  
//  //  // //     // //     // //       //     // //     // 
//  //  // //     // //     // //       //     //        // 
//  //  // //     // ////////  //       //     //  ///////  
//  //  // //     // //   //   //       //     // //        
//  //  // //     // //    //  //       //     // //        
 ///  ///   ///////  //     // //////// ////////  ///////// 
namespace world2 {


/*  Verbatim transcription from Appendix B, Equations of the World Model
    from the book World Dynamics, Second Edition by Jay W. Forrester, 1973.


            *         WORLD DYNAMICS W5
     1      L    P.K=P.J+(DT)(BR.JK-DR.JK)
     1.1    N    P=PI
     1.2    C    PI=1.65E9
     2      R    BR.KL=(P.K)(CLIP(BRN,BRN1,SWT1,TIME.K))(BRFM.K)(BRMM.K)(BRCM.K)(BR
     2.1    X    PM.K)
     2.2    C    BRN=.04
     2.3    C    BRN1=.04
     2.4    C    SWT1=1970
     3      A    BRMM.K=TABHL(BRMMT,MSL.K,0,5,1) 
     3.1    T    BRMMT=1.2/1/.85/.75/.7/.7
     4      A    MSL.K=ECIR.K/(ECIRN)
     4.1    C    ECIRN=1
     5      A    ECIR.K=(CIR.K)(1-CIAF.K)(NREM.K)/(1-CIAFN)
     6      A    NREM.K=TABLE(NREMT,NRFR.K,0,1,.25)
     6.1    T    NREMT=0/.15/.5/.85/1
     7      A    NRFR.K=NR.K/NRI
     8      L    NR.K=NR.J+(DT)(-NRUR.JK)
     8.1    N    NR=NRI
     8.2    C    NRI=900E9
     9      R    NRUR.KL=(P.K)(CLIP(NRUN,NRUN1,SWT2,TIME.K))*(NRMM.K)
     9.1    C    NRUN=1
     9.2    C    NRUN1=1
     9.3    C    SWT2=1970
            NOTE      EQUATION 42 CONNECTS HERE FROM EQ. 4 TO EQ.9
    10      R    DR.KL=(P.K)(CLIP(DRN,DRN1,SWT3,TIME.K))(DRMM.K)(DRPM.K)(DRFM.K)(DR
            X    CM.K)
    10.2    C    DRN=.028
    10.3    C    DRN1=.028
    10.4    C    SWT3=1970
    11      A    DRMM.K=TABHL(DRMMT,MSL.K,0,5,.5)
    11.1    T    DRMMT=3/1.8/1/.8/.7/.6/.53/.5/.5/.5/.5
    12      A    DRPM.K=TABLE(DRPMT,POLR.K,0,60,10)
    12.1    T    DRPMT=.92/1.3/2/3.2/4.8/6.8/9.2
    13      A    DRFM.K=TABHL(DRFMT,FR.K,0,2,.25)
    13.1    T    DRFMT=30/3/2/1.4/1/.7/.6/.5/.5
    14      A    DRCM.K=TABLE(DRCMT,CR.K,0,5,1)
    14.1    T    DRCMT=.9/1/1.2/1.5/1.9/3
    15      A    CR.K=(P.K)/(LA*PDN)
    15.1    C    LA=135E6
    15.2    C    PDN=26.5
    16      A    BRCM.K=TABLE(BRCMT,CR.K,0,5,1)
    16.1    T    BRCMT=1.05/1/.9/.7/.6/.55
    17      A    BRFM.K=TABHL(BRFMT,FR.K,0,4,1)
    17.1    T    BRFMT=0/1/1.6/1.9/2
    18      A    BRPM.K=TABLE(BRPMT,POLR.K,0,60,10)
    18.1    T    BRPMT=1.02/.9/.7/.4/.25/.15/.1
    19      A    FR.K=(FPCI.K)(FCM.K)(FPM.K)(CLIP(FC,FC1,SWT7,TIME.K))/FN
    19.1    C    FC=1
    19.2    C    FC1=1
    19.3    C    FN=1
    19.4    C    SWT7=1970
    20      A    FCM.K=TABLE(FCMT,CR.K,0,5,1)
    20.1    T    FCMT=2.4/1/.6/.4/.3/.2
    21      A    FPCI.K=TABHL(FPCIT,CIRA.K,0,6,1)
    21.1    T    FPCIT=.5/1/1.4/1.7/1.9/2.05/2.2
    22      A    CIRA.K=(CIR.K)(CIAF.K)/CIAFN
    22.1    C    CIAFN=.3
    23      A    CIR.K=CI.K/P.K
    24      L    CI.K=CI.J+(DT)(CIG.JK-CID.JK)
    24.1    N    CI=CII
    24.2    C    CII=.4E9
    25      R    CIG.KL=(P.K)(CIM.K)(CLIP(CIGN,CIGN1,SWT4,TIME.K))
    25.1    C    CIGN=.05
    25.2    C    CIGN1=.05
    25.3    C    SWT4=1970
    26      A    CIM.K=TABHL(CIMT,MSL.K,0,5,1)
    26.1    T    CIMT=.1/1/1.8/2.4/2.8/3
    27      R    CID.KL=(CI.K)(CLIP(CIDN,CIDN1,SWT5,TIME.K))
    27.1    C    CIDN=.025
    27.2    C    CIDN1=.025
    27.3    C    SWT5=1970
    28      A    FPM.K=TABLE(FPMT,POLR.K,0,60,10)
    28.1    T    FPMT=1.02/.9/.65/.35/.2/.1/.05
    29      A    POLR.K=POL.K/POLS
    29.1    C    POLS=3.6E9
    30      L    POL.K=POL.J+(DT)(POLG.JK-POLA.JK)
    30.1    N    POL=POLI
    30.2    C    POLI=.2E9
    31      R    POLG.KL=(P.K)(CLIP(POLN,POLN1,SWT6,TIME.K))(POLCM.K)
    31.1    C    POLN=1
    31.2    C    POLN1=1
    31.3    C    SWT6=1970
    32      A    POLCM.K=TABHL(POLCMT,CIR.K,0,5,1)
    32.1    T    POLCMT=.05/1/3/5.4/7.4/8
    33      R    POLA.KL=POL.K/POLAT.K
    34      A    POLAT.K=TABLE(POLATT,POLR.K,0,60,10)
    34.1    T    POLATT=.6/2.5/5/8/11.5/15.5/20
    35      L    CIAF.K=CIAF.J+(DT/CIAFT)((CFIFR.J*CIQR.J)-CIAF.J)
    35.1    N    CIAF=CIAFI
    35.2    C    CIAFI=.2
    35.3    C    CIAFT=15
    36      A    CFIFR.K=TABHL(CFIFRT,FR.K,0,2,.5)
    36.1    T    CFIFRT=1/.6/.3/.15/.1
    37      A    QL.K=(QLS)(QLM.K)(QLC.K)(QLF.K)(QLP.K)
    37.1    C    QLS=1
    38      A    QLM.K=TABHL(QLMT,MSL.K,0,5,1)
    38.1    T    QLMT=.2/1/1.7/2.3/2.7/2.9
    39      A    QLC.K=TABLE(QLCT,CR.K,0,5,.5)
    39.1    T    QLCT=2/1.3/1/.75/.55/.45/.38/.3/.25/.22/.2
    40      A    QLF.K=TABHL(QLFT,FR.K,0,4,1)
    40.1    T    QLFT=0/1/1.8/2.4/2.7
    41      A    QLP.K=TABLE(QLPT,POLR.K,0,60,10)
    41.1    T    QLPT=1.04/.85/.6/.3/.15/.05/.02
            NOTE      EQUATION 42 LOCATED BETWEEN EQ. 4 AND 9.
    42      A    NRMM.K=TABHL(NRMMT,MSL.K,0,10,1)
    42.1    T    NRMMT=0/1/1.8/2.4/2.9/3.3/3.6/3.8/3.9/3.95/4
            NOTE      INPUT FROM EQN. 38 AND 40 TO EQN. 35
    43      A    CIQR.K=TABHL(CIQRT,QLM.K/QLF.K,0,2,.5)
    43.1    T    CIQRT=.7/.8/1/1.5/2
            NOTE
            NOTE CONTROL CARDS
            NOTE
    43.5    C    DT=.2
    43.6    C    LENGTH=2100
    43.7    N    TIME=1900
    44      A    PRTPER.K=CLIP(PRTP1,PRTP2,PRSWT,TIME.K)
    44.1    C    PRTP1=0
    44.2    C    PRTP2=0
    44.3    C    PRSWT=0
    45      A    PLTPER.K=CLIP(PLTP1,PLTP2,PLSWT,TIME.K)
    45.1    C    PLTP1=4
    45.2    C    PLTP2=4
    45.3    C    PLSWT=0
            PLOT P=P(0,8E9)/POLR=2(0,40)/CI=C(0,20E9)/QL=Q(0,2)/NR=N(0,1000E9)
            NOTE PLOT FR=F,MSL=M,QLC=4,QLP=5(0,2)/CIAF=A(.2,.6) 
            RUN  ORIG


    We're going to re-imagine this DYNAMO code in C++. Some excerpts from the
    Dynamo User's Manual [by Alexander L. Pugh, III, Second Edition, 1963]
    will help us understand the meaning of the above code:


    1.1 Explanation of Time Notation and Variable Types
    ---------------------------------------------------
    The time notation and the types of variables are discussed in detail in
    Chapter 7 of Industrial Dynamics but will be described briefly here.

    The basis for the time notation is the procedure by which the computer
    calculates the results, which is to move through TIME in discrete steps
    and calculate all the variables at each step. Figure 1-1 shows the
    procedure graphically.


          |              |              |
          |<---- JK ---->|<---- KL ---->|
          |   interval   |   interval   |
          |              |              |
        --+--------------+--------------+------->
          |              |              |    TIME
          J              K              L
                      (Present
                        TIME)

               Figure 1-1 Time Notation.


    The TIME for which the calculations are currently being made is called
    TIME K. The previous TIME for which calculations were made is called J,
    and the next instant for which calculations will be made is L. The
    intervals between these times are called JK and KL. The length of these
    intervals is called DT.

    The names of instants and the intervals are used to specify when a
    double is calculated and when the quantities used in the calculation
    were previously calculated. Once all the variables have been calculated
    for the instant K and the interval KL, the computer moves forward one
    time step, and the values that were associated with TIME K are now
    related to TIME J.
    -- page 4

    ...

    There are three principal types of variables in DYNAMO: levels, rates,
    and auxiliaries.

    Level
    -----
    A level, which is calculated at TIME K, is a double that depends upon
    its value at TIME J and on other quantities at that TIME or in the JK
    interval. Inventory is an example in that the inventory today is equal
    to the inventory at an earlier time plus what has been received minus
    what has been shipped during the interim.

    Rate
    ----
    The decisions in business and economic models are rates. Rates are the
    flows of tangible things from one level to the next. They are computed
    at TIME K for the interval KL from levels and auxiliaries at TIME K and
    occasionally from rates in the JK interval.

    There is some confusion about rates in that all quantities with the
    dimension of something per unit time are not necessarily rates. In
    particular, averaged or smoothed rates, which are measured in something
    per unit time, are levels.

    Auxiliary
    ---------
    Auxiliaries are variables that are introduced to simplify the algebraic
    complexity of rate equations. They generally have some physical meaning
    and consequently simplify the understanding of the model. They are
    computed at TIME K from levels and other auxiliaries at the same time
    and occasionally from rates in the JK interval. By their nature they can
    be eliminated by substitution into the rate equations.

    Order of Computation
    --------------------
    The order of computation at TIME K is first levels, which are based on
    quantities from TIME J and the JK interval, next auxiliaries, which are
    based on levels and auxiliaries computed earlier at TIME K and on rates
    JK, and finally rates, which are based on levels and auxiliaries from K
    and other rates from JK.
    -- page 5

    ...

    3.3 Running Phase
    -----------------
    After all the instructions for a model have been generated, it is ready
    to be run. The first step is to initialize those quantities that require
    initial values. Thus, DYNAMO

        1. Sets TIME equal to 0
        2. Temporarily sets all computed constants, levels, auxiliaries,
            rates, and supplementaries equal to 0
        3. Loads all the boxcar trains with the initial values provided by
            the special C cards
        4. Computes the computed constants, and the initial values of levels
            and of those auxiliaries and rates for which initial values were
            provided (that is, computes all N-type equations)

    Now the model is initialized and ready to generate data. As the initial
    values of the levels have just been computed, DYNAMO starts with the
    auxiliaries or step 4 below. A complete time step consists of the
    following operations:

        1. The differences between levels at time .J and time .K are computed.
        2. The differences are added to the levels at time .J to form the
            levels at time .K.
        3. DT is added to TIME.J to form TIME.K.
        4. Auxiliaries are computed for time .K.
        5. Rates for the interval .KL are computed and stored in the extra
            location for rates. Thus the rates for the interval .JK are
            undisturbed and available for the calculation of other rates.
        6. TIME is tested to determine whether this is an output time.
            (TIME=0 is, of course, an output time.) If so, the following
            three substeps are executed. If not, only substep b is executed,
            and DYNAMO proceeds on to step 7.

            a. The supplementaries are computed using the .JK values of the
                rates.
            b. The rates just computed for the .KL interval are shifted from
                their extra location to the normal location for rates.
            c. The present values of those quantities that are to be printed
                or plotted are saved on magnetic tape for later processing.

    In order that DYNAMO can choose a scaling factor for printing and the
    upper and lower scale for plotting, DYNAMO must have the maximum and
    minimum values of all the quantities being saved on magnetic tape.
    Therefore, as each double is located and saved, it is compared with
    the past maximum and minimum for that double. If either the past
    maximum or minimum is exceeded, that value is replaced with the current
    value of the double.

        7. TIME is tested to determine whether the run is completed. If so,
            DYNAMO proceeds to the printing phase. If not, DYNAMO goes on
            to the next step.
        8. Each boxcar train is tested and shifted if this TIME is
            appropriate. (All trains are shifted at TIME=0.)

    At this point, time .L becomes time .K, time .K becomes time .J, and
    DYNAMO starts the time-step procedure again with step 1 just mentioned.
    -- page 48

*/




class world {

    // numbers in square brackets refer to the line numbers of
    // Forrester's original World2 DYNAMO code

public:
    struct constants {
        double brn      = .04;      //[2.2]     birth rate normal (fraction/year)
        double brn1     = .04;      //[2.3]     birth rate normal no. 1 (fraction/year)
        double ciafi    = .2;       //[35.2]    capital-investment-in-agriculture-fraction initial (dimensionless)
        double ciafn    = .3;       //[22.1]    capital-investment-in-agriculture fraction normal (dimensionless)
        double ciaft    = 15;       //[35.3]    capital-investment-in-agriculture-fraction adjustment time (years)
        double cidn     = .025;     //[27.1]    capital-investment discard normal (fraction/year)
        double cidn1    = .025;     //[27.2]    capital-investment discard normal no. 1 (fraction/year)
        double cign     = .05;      //[25.1]    capital-investment generation normal (capital units/person/year)
        double cign1    = .05;      //[25.2]    capital-investment generation normal no. 1 (capital units/person/year)
        double cii      = .4E9;     //[24.2]    capital-investment, initial (capital units)
        double drn      = .028;     //[10.2]    death rate normal (fraction/year)
        double drn1     = .028;     //[10.3]    death rate normal no. 1 (fraction/year)
        double ecirn    = 1;        //[4.1]     effective-capital-investment ratio normal (capital units/person)
        double fc       = 1;        //[19.1]    food coefficient (dimensionless)
        double fc1      = 1;        //[19.2]    food coefficient no. 1 (dimensionless)
        double fn       = 1;        //[19.3]    food normal (food units/person/year)
        double la       = 135E6;    //[15.1]    land area (square kilometers)
        double nri      = 900E9;    //[8.2]     natural resources, initial (natural resource units)
        double nrun     = 1;        //[9.1]     natural-resource usage normal (natural resource units/person/year)
        double nrun1    = 1;        //[9.2]     natural-resource usage normal no. 1 (natural resource units/person/year)
        double pdn      = 26.5;     //[15.2]    population density normal (people/square kilometer)
        double pi       = 1.65E9;   //[1.1]     population, initial (people)
        double poli     = .2E9;     //[30.2]    pollution, initial (pollution units)
        double poln     = 1;        //[31.1]    pollution normal (pollution units/person/year)
        double poln1    = 1;        //[31.2]    pollution normal no. 1 (pollution units/person/year)
        double pols     = 3.6E9;    //[29.1]    pollution standard (pollution units)
        double qls      = 1;        //[37.1]    quality-of-life standard (satisfaction units)
        double swt1     = 1970;     //[2.4]     switch time no. 1 for brn (years)
        double swt2     = 1970;     //[9.3]     switch time no. 2 for nrun (years)
        double swt3     = 1970;     //[10.4]    switch time no. 3 for drn (years)
        double swt4     = 1970;     //[25.3]    switch time no. 4 for cign (years)
        double swt5     = 1970;     //[27.3]    switch time no. 5 for cidn (years)
        double swt6     = 1970;     //[31.3]    switch time no. 6 for poln (years)
        double swt7     = 1970;     //[19.4]    switch time no. 7 for fc (years)

        double time     = 1900;     //[43.7]    calendar time (years)
        double dt       = 0.2;      //[43.5]    delta time (years)
        double endtime  = 2100;     // when time has this value the run should terminate
    };

    struct variables {
        // levels
        double ci   = 0;    // capital-investment (capital units)
        double ciaf = 0;    // capital-investment-in-agriculture fraction
        double nr   = 0;    // natural resources (natural resource units)
        double p    = 0;    // population
        double pol  = 0;    // pollution (pollution units)

        // rates
        double br   = 0;    // birth rate (people/year)
        double cid  = 0;    // capital-investment discard (capital units/year)
        double cig  = 0;    // capital-investment generation (capital units/year)
        double dr   = 0;    // death rate (people/year)
        double nrur = 0;    // natural-resource-usage rate (natural resource units/year)
        double pola = 0;    // pollution absorption (pollution units/year)
        double polg = 0;    // pollution generation (pollution units/year)

        // auxilaries
        double brcm = 0;    // birth-rate-from-crowding multiplier
        double brfm = 0;    // birth-rate-from-food multiplier
        double brmm = 0;    // birth-rate-from-material multiplier
        double brpm = 0;    // birth-rate-from-pollution multiplier
        double cfifr = 0;   // capital fraction indicated by food ratio
        double cim  = 0;    // capital-investment multiplier
        double ciqr = 0;    // capital-investment-from-quality ratio
        double cir  = 0;    // capital-investment ratio (capital units/person)
        double cira = 0;    // capital-investment ratio in agriculture (capital units/person)
        double cr   = 0;    // crowding ratio
        double drcm = 0;    // death-rate-from-crowding multiplier
        double drfm = 0;    // death-rate-from-food multiplier
        double drmm = 0;    // death-rate-from-material multiplier
        double drpm = 0;    // death-rate-from-pollution multiplier
        double ecir = 0;    // effective-capital-investment ratio (capital units/person)
        double fcm  = 0;    // food-from-crowding multiplier
        double fpci = 0;    // food potential from capital investment (food units/person/year)
        double fpm  = 0;    // food-from-pollution multiplier
        double fr   = 0;    // food ratio
        double msl  = 0;    // material standard of living
        double nrem = 0;    // natural-resource-extraction multiplier
        double nrfr = 0;    // natural-resource fraction remaining
        double nrmm = 0;    // natural-resource-from-material multiplier
        double polat = 0;   // pollution-absorption time (years)
        double polcm = 0;   // pollution-from-capital multiplier
        double polr = 0;    // pollution ratio
        double ql   = 0;    // quality of life
        double qlc  = 0;    // quality of life from crowding
        double qlf  = 0;    // quality of life from food
        double qlm  = 0;    // quality of life from material
        double qlp  = 0;    // quality of life from pollution

        double time = 0;    // calendar time (years)
    };


    world(const constants & c)
        : c(c), time_j_exists_(false)
    {
    }

    bool run_complete() const
    {
        return time_j_exists_ && j.time > c.endtime;
    }

    // return a reference to variables calculated for time .K
    const variables & tick()
    {
        variables k;

        if (time_j_exists_) {
            // calculate levels at time .K (note that .JK rates are shown here as j.xxx)
            k.p     = j.p + c.dt * (j.br - j.dr);       //[1]
            k.nr    = j.nr + c.dt * -j.nrur;            //[8]
            k.ci    = j.ci + c.dt * (j.cig - j.cid);    //[24]
            k.pol   = j.pol + c.dt * (j.polg - j.pola); //[30]
            k.ciaf  = j.ciaf + (c.dt / c.ciaft) * ((j.cfifr * j.ciqr) - j.ciaf); //[35]

            k.time  = j.time + c.dt;
        }
        else {
            // set levels to initial state
            k.p     = c.pi;
            k.nr    = c.nri;
            k.ci    = c.cii;
            k.pol   = c.poli;
            k.ciaf  = c.ciafi;

            k.time  = c.time;
        }

        // compute auxiliaries for time .K (reordered for dependencies)
        k.nrfr  = k.nr / c.nri; //[7]
        k.nrem  = dynamo::table({ 0, .15, .5, .85, 1 }, k.nrfr, 0, 1, .25); //[6, 6.1]
        k.cir   = k.ci / k.p; //[23]
        k.ecir  = k.cir * (1 - k.ciaf) * k.nrem / (1 - c.ciafn); //[5]
        k.msl   = k.ecir / c.ecirn; //[4]
        k.brmm  = dynamo::tabhl({ 1.2, 1, .85, .75, .7, .7 }, k.msl, 0, 5, 1); //[3, 3.1]
        k.drmm  = dynamo::tabhl({ 3, 1.8, 1, .8, .7, .6, .53, .5, .5, .5, .5 }, k.msl, 0, 5, .5); //[11, 11.1]
        k.cr    = k.p / (c.la * c.pdn); //[15]
        k.drcm  = dynamo::table({ .9, 1, 1.2, 1.5, 1.9, 3 }, k.cr, 0, 5, 1); //[14, 14.1]
        k.brcm  = dynamo::table({ 1.05, 1, .9, .7, .6, .55 }, k.cr, 0, 5, 1); //[16, 16.1]
        k.fcm   = dynamo::table({ 2.4, 1, .6, .4, .3, .2 }, k.cr, 0, 5, 1); //[20, 20.1]
        k.qlc   = dynamo::table({ 2, 1.3, 1, .75, .55, .45, .38, .3, .25, .22, .2 }, k.cr, 0, 5, .5); //[39, 39.1]
        k.cim   = dynamo::tabhl({ .1, 1, 1.8, 2.4, 2.8, 3 }, k.msl, 0, 5, 1); //[26, 26.1]
        k.polr  = k.pol / c.pols; //[29, 29.1]
        k.fpm   = dynamo::table({ 1.02, .9, .65, .35, .2, .1, .05 }, k.polr, 0, 60, 10); //[28, 28.1]
        k.drpm  = dynamo::table({ .92, 1.3, 2, 3.2, 4.8, 6.8, 9.2 }, k.polr, 0, 60, 10); //[12, 12.1]
        k.brpm  = dynamo::table({ 1.02, .9, .7, .4, .25, .15, .1 }, k.polr, 0, 60, 10); //[18, 18.1]
        k.polcm = dynamo::tabhl({ .05, 1, 3, 5.4, 7.4, 8 }, k.cir, 0, 5, 1); //[32, 32.1]
        k.polat = dynamo::table({ .6, 2.5, 5, 8, 11.5, 15.5, 20 }, k.polr, 0, 60, 10); //[34, 34.1]
        k.qlm   = dynamo::tabhl({ .2, 1, 1.7, 2.3, 2.7, 2.9 }, k.msl, 0, 5, 1); //[38, 38.1]
        k.qlp   = dynamo::table({ 1.04, .85, .6, .3, .15, .05, .02 }, k.polr, 0, 60, 10); //[41, 41.1]
        k.nrmm  = dynamo::tabhl({ 0, 1, 1.8, 2.4, 2.9, 3.3, 3.6, 3.8, 3.9, 3.95, 4 }, k.msl, 0, 10, 1); //[42, 42.1]
        k.cira  = k.cir * k.ciaf / c.ciafn; //[22]
        k.fpci  = dynamo::tabhl({ .5, 1, 1.4, 1.7, 1.9, 2.05, 2.2 }, k.cira, 0, 6, 1); //[21, 21.1]
        k.fr    = k.fpci * k.fcm * k.fpm * dynamo::clip(c.fc, c.fc1, c.swt7, k.time) / c.fn; //[19]
        k.drfm  = dynamo::tabhl({ 30, 3, 2, 1.4, 1, .7, .6, .5, .5 }, k.fr, 0, 2, .25); //[13, 13.1]
        k.brfm  = dynamo::tabhl({ 0, 1, 1.6, 1.9, 2 }, k.fr, 0, 4, 1); //[17, 17.1]
        k.cfifr = dynamo::tabhl({ 1, .6, .3, .15, .1 }, k.fr, 0, 2, .5); //[36, 36.1]
        k.qlf   = dynamo::tabhl({ 0, 1, 1.8, 2.4, 2.7 }, k.fr, 0, 4, 1); //[40, 40.1]
        k.ciqr  = dynamo::tabhl({ .7, .8, 1, 1.5, 2 }, k.qlm / k.qlf, 0, 2, .5); //[43, 43.1]
        k.ql    = c.qls * k.qlm * k.qlc * k.qlf * k.qlp; //[37]

        // calculate rates for period .KL (write direct to .JK as no references to .JK are made)
        k.br    = k.p * dynamo::clip(c.brn, c.brn1, c.swt1, k.time) * k.brfm * k.brmm * k.brcm * k.brpm; //[2, 2.1]
        k.nrur  = k.p * dynamo::clip(c.nrun, c.nrun1, c.swt2, k.time) * k.nrmm; //[9]
        k.dr    = k.p * dynamo::clip(c.drn, c.drn1, c.swt3, k.time) * k.drmm * k.drpm * k.drfm * k.drcm; //[10, 10.1]
        k.cig   = k.p * k.cim * dynamo::clip(c.cign, c.cign1, c.swt4, k.time); //[25]
        k.cid   = k.ci * dynamo::clip(c.cidn, c.cidn1, c.swt5, k.time); //[27]
        k.polg  = k.p * dynamo::clip(c.poln, c.poln1, c.swt6, k.time) * k.polcm; //[31]
        k.pola  = k.pol / k.polat; //[33]

        // shift .K to .J for next call to tick()
        j = k;
        time_j_exists_ = true;

        return j; // which on this tick is .K
    }

private:
    constants c;
    variables j;
    bool time_j_exists_ = false;
};



}//namespace world2






 //////   ////////     ///    ////////  //     // 
//    //  //     //   // //   //     // //     // 
//        //     //  //   //  //     // //     // 
//   //// ////////  //     // ////////  ///////// 
//    //  //   //   ///////// //        //     // 
//    //  //    //  //     // //        //     // 
 //////   //     // //     // //        //     // 

using namespace world2;

// Approximate the DYNAMO graphs as shown in Forrester's book.
// This is not a full DYNAMO graph implementation, but is sufficient
// to draw the graphs I want to show here.
class graph {
public:
    graph(const world::constants & c)
        : w_(c)
    {}


    void plot(
        double world::variables::* vptr,
        const char * name,
        const char symbol,
        double low,
        double high)
    {
        plotvars_.push_back({ vptr, symbol, low, high });

        if (!ledgend_.empty())
            ledgend_ += ',';
        ledgend_ += name;
        ledgend_ += '=';
        ledgend_ += symbol;
    }

    std::string run()
    {
        const size_t graph_width = default_graph_width;
        const size_t x_label_width = 8;
        const size_t dash_interval = 10;
        const size_t dot_interval = 15;
        size_t row_counter = 0;

        std::vector<std::string> lines;

        for (size_t tick = 0; !w_.run_complete(); ++tick) {
            const world::variables & vars = w_.tick();
            if (tick % 20 == 0) {
                std::string line(graph_width+1, ' ');
                std::string x_label(x_label_width, ' ');

                for (size_t i = 0; i < graph_width+1; i += dot_interval)
                    line[i] = '.';
                if (row_counter++ % dash_interval == 0) {
                    for (size_t i = 0; i < graph_width+1; i += 2) {
                        line[i] = '-';
                        if (i < graph_width)
                            line[i + 1] = ' ';
                    }
                    char buf[100];
                    snprintf(buf, sizeof(buf), "%7.0f.", vars.time);
                    x_label = buf;
                }

                std::map<char, std::string> intersects;

                for (const plotvar & pv : plotvars_) {
                    const int y = calc_y(vars.*(pv.vptr), pv.low, pv.high, graph_width);
                    if (y < graph_width+1) {
                        if (line[y] == ' ' || line[y] == '-' || line[y] == '.')
                            line[y] = pv.symbol;
                        else {
                            intersects[(line[y])] += pv.symbol;
                        }
                    }
                }

                if (!intersects.empty()) {
                    line += ' ';
                    bool need_comma = false;
                    for (auto & i : intersects) {
                        if (need_comma)
                            line += ',';
                        line += i.first;
                        line += i.second;
                        need_comma = true;
                    }
                }

                lines.push_back(x_label + line);
            }
        }

        std::vector<std::string> y_scale;
        for (const plotvar & pv : plotvars_) {
            std::string scale(x_label_width - 2, ' ');
            scale += pv.symbol;
            scale += ' ';
            const size_t steps = static_cast<size_t>(graph_width / dot_interval);
            const double step = (pv.high - pv.low) / steps;
            double label = pv.low;
            for (size_t i = 0; i < steps; ++i) {
                char buf[100];
                snprintf(buf, sizeof(buf), "%-*s",
                    static_cast<int>(dot_interval),
                    numeric_fmt(label).c_str());
                scale += buf;
                label += step;
            }
            scale += numeric_fmt(pv.high);
            y_scale.push_back(scale);
        }

        return ledgend_ + "\n\n" + join(y_scale, "\n") + "\n" + join(lines, "\n");
    }

private:
    const static size_t default_graph_width = 60;

    world w_;

    struct plotvar {
        double world::variables:: * vptr;   // value to plot
        const char symbol;                  // symbol used to represent value
        double low, high;                   // y-axis lower and upper bounds
    };
    std::vector<plotvar> plotvars_;
    std::string ledgend_;

    static std::string join(const std::vector<std::string> & strings, const std::string & joiner)
    {
        std::string result;
        for (auto & s : strings) {
            if (!result.empty())
                result += joiner;
            result += s;
        }
        return result;
    }

    static int calc_y(double value, double scale_lo, double scale_hi, int num_divisions)
    {
        return static_cast<int>(std::round((value - scale_lo) / (scale_hi - scale_lo) * num_divisions));
    }

    static std::string numeric_fmt(double d)
    {
        char buf[100] = "";
        if (std::fabs(d) > 1e12)
            snprintf(buf, sizeof(buf), "%g", d);
        else {
            char magnitude = '\0';
            if (std::fabs(d) >= 1e9) {
                d /= 1e9;
                magnitude = 'B';
            }
            else if (std::fabs(d) >= 1e6) {
                d /= 1e6;
                magnitude = 'M';
            }
            else if (std::fabs(d) >= 1e3) {
                d /= 1e3;
                magnitude = 'T';
            }
            snprintf(buf, sizeof(buf), "%.3f", d);

            char * cp = buf + std::strlen(buf) - 1;
            while (*cp == '0')
                *cp-- = '\0';
            *++cp = magnitude;
            *++cp = '\0';
        }
        return buf;
    }

    friend void test();
};




void fig_41()
{
    graph g({});
    //PLOT P=P(0,8E9)/POLR=2(0,40)/CI=C(0,20E9)/QL=Q(0,2)/NR=N(0,1000E9)
    g.plot(&world::variables::p,    "P",    'P', 0, 8E9);
    g.plot(&world::variables::polr, "POLR", '2', 0, 40);
    g.plot(&world::variables::ci,   "CI",   'C', 0, 20E9);
    g.plot(&world::variables::ql,   "QL",   'Q', 0, 2);
    g.plot(&world::variables::nr,   "NR",   'N', 0, 1000E9);
    std::cout
        << "\n\n"
        << g.run()
        << "\n\n"
        << "    World Dynamics, Figure 4-1, ORIG-A\n"
        << "    \"Basic behavior of the world model, showing the mode in which\n"
        << "     industrialization and population are suppressed by falling\n"
        << "     natural resources.\"\n"
        << "    [P=Population, 2=Pollution, C=Capital Investment, Q=Quality of Life,\n"
        << "     N=Natural Resources]\n";
}

void fig_42()
{
    graph g({});
    //PLOT FR=F,MSL=M,QLC=4,QLP=5(0,2)/CIAF=A(.2,.6) 
    g.plot(&world::variables::fr,   "FR",   'F', 0, 2);
    g.plot(&world::variables::msl,  "MSL",  'M', 0, 2);
    g.plot(&world::variables::qlc,  "QLC",  '4', 0, 2);
    g.plot(&world::variables::qlp,  "QLP",  '5', 0, 2);
    g.plot(&world::variables::ciaf, "CIAF", 'A', .2, .6);
    std::cout
        << "\n\n"
        << g.run()
        << "\n\n"
        << "    World Dynamics, Figure 4-2, ORIG-B\n"
        << "    \"Original model as in Figure 4-1. Material standard of living reaches\n"
        << "     a maximum and then declines as natural resources are depleted.\"\n"
        << "    [F=Food ratio, M=Material standard of living, 4=Quality of Life from\n"
        << "     crowding, 5=Quality of life from pollution, A=Capital-investment-in-\n"
        << "     agriculture fraction]\n";
}

void fig_43()
{
    graph g({});
    g.plot(&world::variables::nr,   "NR",   'N', 0, 1e12);
    g.plot(&world::variables::nrur, "NRUR", 'U', 0, 8e9);
    std::cout
        << "\n\n"
        << g.run()
        << "\n\n"
        << "    World Dynamics, Figure 4-3, ORIG-C\n"
        << "    \"Original model as in Figure 4-1. Natural-resource-usage rate reaches\n"
        << "     a peak about year 2010 and declines as natural resources, population,\n"
        << "     and capital investment decline.\"\n"
        << "    [N=Natural resources, U=Natural-resource-usage rate]\n";
}

void fig_44()
{
    graph g({});
    g.plot(&world::variables::ci,   "CI",   'C', 0, 20E9);
    g.plot(&world::variables::cig,  "CIG",  'G', 0, 400E6);
    g.plot(&world::variables::cid,  "CID",  'D', 0, 400E6);
    std::cout
        << "\n\n"
        << g.run()
        << "\n\n"
        << "    World Dynamics, Figure 4-4, ORIG-D\n"
        << "    \"Original model as in Figure 4-1. The rate of capital-investment generation\n"
        << "     declines after 2010 but does not fall below the rate of capital-investment\n"
        << "     discard until 2040, at which time the level of capital investment begins\n"
        << "     to decline.\"\n"
        << "    [C=Capital investment, G=Capital-investment generation, D=Capital-investment\n"
        << "     discard]\n";
}

void fig_45()
{
    world::constants c;
    c.nrun1 = 0.25;
    graph g(c);
    g.plot(&world::variables::p,    "P",    'P', 0, 8E9);
    g.plot(&world::variables::polr, "POLR", '2', 0, 40);
    g.plot(&world::variables::ci,   "CI",   'C', 0, 20E9);
    g.plot(&world::variables::ql,   "QL",   'Q', 0, 2);
    g.plot(&world::variables::nr,   "NR",   'N', 0, 1000E9);
    std::cout
        << "\n\n"
        << g.run()
        << "\n\n"
        << "    World Dynamics, Figure 4-5, Original NRUN1=1.0, present NRUN1=0.25\n"
        << "    \"Reduced usage rate of natural resources leads to a pollution crisis.\"\n"
        << "    [P=Population, 2=Pollution, C=Capital Investment, Q=Quality of Life,\n"
        << "     N=Natural Resources]\n";
}

void fig_46()
{
    world::constants c;
    c.nrun1 = 0.25;
    graph g(c);
    g.plot(&world::variables::fr,   "FR",   'F', 0, 2);
    g.plot(&world::variables::msl,  "MSL",  'M', 0, 2);
    g.plot(&world::variables::qlc,  "QLC",  '4', 0, 2);
    g.plot(&world::variables::qlp,  "QLP",  '5', 0, 2);
    g.plot(&world::variables::ciaf, "CIAF", 'A', .2, .6);
    std::cout
        << "\n\n"
        << g.run()
        << "\n\n"
        << "    World Dynamics, Figure 4-6, Original NRUN1=1.0, present NRUN1=0.25\n"
        << "    \"System ratios during the pollution mode of growth suppression.\"\n"
        << "    [F=Food ratio, M=Material standard of living, 4=Quality of Life from\n"
        << "     crowding, 5=Quality of life from pollution, A=Capital-investment-in-\n"
        << "     agriculture fraction]\n";
}

void fig_47()
{
    world::constants c;
    c.nrun1 = 0.25;
    graph g(c);
    g.plot(&world::variables::polr, "POLR", '2', 0, 40);
    g.plot(&world::variables::polat,"POLAT",'T', 0, 16);
    g.plot(&world::variables::polg, "POLG", 'G', 0, 20E9);
    g.plot(&world::variables::pola, "POLA", 'A', 0, 20E9);
    std::cout
        << "\n\n"
        << g.run()
        << "\n\n"
        << "    World Dynamics, Figure 4-7, Original NRUN1=1.0, present NRUN1=0.25\n"
        << "    \"Dynamics of the pollution sector. A positive-feedback growth\n"
        << "     in pollution occurs when the pollution-absorption time increases\n"
        << "     faster than the pollution.\"\n"
        << "    [2=Pollution, T=Pollution-absorption time, G=Pollution generation,\n"
        << "     A=Pollution absorption]\n";
}

// there are more graphs in Forrester's book, but I'm not going
// to recreate them all here



void test()
{
    const std::vector<double> t1{ 1.0, 2.0 };
    TEST_EQUAL_DOUBLE(dynamo::table(t1, 3.0, 3.0, 4.0, 1.0), 1.0);
    TEST_EQUAL_DOUBLE(dynamo::table(t1, 4.0, 3.0, 4.0, 1.0), 2.0);
    TEST_EQUAL_DOUBLE(dynamo::table(t1, 3.5, 3.0, 4.0, 1.0), 1.5);
    TEST_EQUAL_DOUBLE(dynamo::table(t1, 3.75, 3.0, 4.0, 1.0), 1.75);

    TEST_EQUAL_DOUBLE(dynamo::table(t1, 3.0, 4.0, 3.0, -1.0), 2.0);
    TEST_EQUAL_DOUBLE(dynamo::table(t1, 4.0, 4.0, 3.0, -1.0), 1.0);
    TEST_EQUAL_DOUBLE(dynamo::table(t1, 3.5, 4.0, 3.0, -1.0), 1.5);
    TEST_EQUAL_DOUBLE(dynamo::table(t1, 3.75, 4.0, 3.0, -1.0), 1.25);

    TEST_EQUAL_DOUBLE(dynamo::table(t1, -0.5, -0.5, 0.5, 1.0), 1.0);
    TEST_EQUAL_DOUBLE(dynamo::table(t1, 0.5, -0.5, 0.5, 1.0), 2.0);
    TEST_EQUAL_DOUBLE(dynamo::table(t1, 0, -0.5, 0.5, 1.0), 1.5);
    TEST_EQUAL_DOUBLE(dynamo::table(t1, 0.25, -0.5, 0.5, 1.0), 1.75);

    TEST_EQUAL_DOUBLE(dynamo::table(t1, -0.5, 0.5, -0.5, -1.0), 2.0);
    TEST_EQUAL_DOUBLE(dynamo::table(t1, 0.5, 0.5, -0.5, -1.0), 1.0);
    TEST_EQUAL_DOUBLE(dynamo::table(t1, 0, 0.5, -0.5, -1.0), 1.5);
    TEST_EQUAL_DOUBLE(dynamo::table(t1, 0.25, 0.5, -0.5, -1.0), 1.25);

    const std::vector<double> t2{ -1.0, 1.0 };
    TEST_EQUAL_DOUBLE(dynamo::table(t2, 3.0, 3.0, 4.0, 1.0), -1.0);
    TEST_EQUAL_DOUBLE(dynamo::table(t2, 4.0, 3.0, 4.0, 1.0), 1.0);
    TEST_EQUAL_DOUBLE(dynamo::table(t2, 3.5, 3.0, 4.0, 1.0), 0);
    TEST_EQUAL_DOUBLE(dynamo::table(t2, 3.75, 3.0, 4.0, 1.0), .5);

    TEST_EQUAL_DOUBLE(dynamo::table(t2, 3.0, 4.0, 3.0, -1.0), 1.0);
    TEST_EQUAL_DOUBLE(dynamo::table(t2, 4.0, 4.0, 3.0, -1.0), -1.0);
    TEST_EQUAL_DOUBLE(dynamo::table(t2, 3.5, 4.0, 3.0, -1.0), 0);
    TEST_EQUAL_DOUBLE(dynamo::table(t2, 3.75, 4.0, 3.0, -1.0), -0.5);

    TEST_EQUAL_DOUBLE(dynamo::table({0.479425, 0.84147099}, 0.632, 0.5, 1.0, 0.5), 0.57500514136);

    const std::vector<double> t3{ 1.04, .85, .6, .3, .15, .05, .02 };
    TEST_EQUAL_DOUBLE(dynamo::table(t3, 0, 0, 60, 10), 1.04);
    TEST_EQUAL_DOUBLE(dynamo::table(t3, 1, 0, 60, 10), 1.0210000000000001);
    TEST_EQUAL_DOUBLE(dynamo::table(t3, 10, 0, 60, 10), .85);
    TEST_EQUAL_DOUBLE(dynamo::table(t3, 20, 0, 60, 10), .6);
    TEST_EQUAL_DOUBLE(dynamo::table(t3, 25, 0, 60, 10), .45);
    TEST_EQUAL_DOUBLE(dynamo::table(t3, 30, 0, 60, 10), .3);
    TEST_EQUAL_DOUBLE(dynamo::table(t3, 40, 0, 60, 10), .15);
    TEST_EQUAL_DOUBLE(dynamo::table(t3, 50, 0, 60, 10), .05);
    TEST_EQUAL_DOUBLE(dynamo::table(t3, 59, 0, 60, 10), .023);
    TEST_EQUAL_DOUBLE(dynamo::table(t3, 60, 0, 60, 10), .02);

    // test the graphs are identical to those printed in the World Dynamics book
    {
        // (actually, just check two curves in one graph)

        //  character positions eyeballed from Figure 4-1 in the World Dynamics book
        const int expected_p_y[51] = {              // scale 0 .. 8e9
            12, 12, 12, 13, 13, 14, 15, 16, 17, 18,
            19, 20, 21, 22, 23, 25, 26, 27, 28, 30,
            31, 32, 34, 35, 36, 37, 38, 39, 39, 40,
            40, 40, 39, 39, 39, 38, 37, 37, 36, 35,
            34, 34, 33, 32, 32, 31, 30, 30, 29, 28,
            28
        };
        const int expected_polr_y[51] = {           // scale 0 .. 40
             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
             0,  1,  1,  1,  1,  1,  1,  1,  1,  2,
             2,  2,  2,  3,  3,  4,  4,  4,  5,  5,
             6,  6,  7,  7,  8,  8,  8,  9,  9,  8,
             8,  8,  8,  7,  7,  6,  6,  5,  5,  4,
             4
        };

        world w({});
        size_t x = 0;
        for (size_t tick = 0; !w.run_complete(); ++tick) {
            const world::variables & vars = w.tick();
            if (tick % 20 == 0) {
                const int p_y = graph::calc_y(vars.p, 0, 8e9, graph::default_graph_width);
                TEST_EQUAL(p_y, expected_p_y[x]);
                const int polr_y = graph::calc_y(vars.polr, 0, 40, graph::default_graph_width);
                TEST_EQUAL(polr_y, expected_polr_y[x]);
                ++x;
            }
        }
   }

    TEST_EQUAL(graph::numeric_fmt(0.0),       "0.");
    TEST_EQUAL(graph::numeric_fmt(1.0/3.0),   "0.333");
    TEST_EQUAL(graph::numeric_fmt(0.5),       "0.5");
    TEST_EQUAL(graph::numeric_fmt(1.0),       "1.");
    TEST_EQUAL(graph::numeric_fmt(2.0),       "2.");
    TEST_EQUAL(graph::numeric_fmt(5.0),       "5.");
    TEST_EQUAL(graph::numeric_fmt(20.0),      "20.");
    TEST_EQUAL(graph::numeric_fmt(250.0),     "250.");
    TEST_EQUAL(graph::numeric_fmt(200e6),     "200.M");
    TEST_EQUAL(graph::numeric_fmt(10000e6/3), "3.333B");
    TEST_EQUAL(graph::numeric_fmt(10000e6),   "10.B");
    TEST_EQUAL(graph::numeric_fmt(250e9),     "250.B");
    TEST_EQUAL(graph::numeric_fmt(1000e9),    "1000.B");
}




int main(int, const char * [])
{
    try {
        test();

        fig_41();
        fig_42();
        fig_43();
        fig_44();
        fig_45();
        fig_46();
        fig_47();
    }
    catch (const std::exception & e) {
        std::cerr << "exception: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    catch (...) {
        std::cerr << "exception" << std::endl;
        return EXIT_FAILURE;
    }
}
