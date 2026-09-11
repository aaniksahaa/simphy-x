#include "sx_coal.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sys.h>
#include <mpfr.h>

double sx_logscale_fact(int n)
{
    double x = 0;
    static const double lf[255] = {
        0.000000000000000,
        0.000000000000000,
        0.693147180559945,
        1.791759469228055,
        3.178053830347946,
        4.787491742782046,
        6.579251212010101,
        8.525161361065415,
        10.604602902745251,
        12.801827480081469,
        15.104412573075516,
        17.502307845873887,
        19.987214495661885,
        22.552163853123421,
        25.191221182738683,
        27.899271383840894,
        30.671860106080675,
        33.505073450136891,
        36.395445208033053,
        39.339884187199495,
        42.335616460753485,
        45.380138898476908,
        48.471181351835227,
        51.606675567764377,
        54.784729398112319,
        58.003605222980518,
        61.261701761002001,
        64.557538627006323,
        67.889743137181526,
        71.257038967168000,
        74.658236348830158,
        78.092223553315307,
        81.557959456115029,
        85.054467017581516,
        88.580827542197682,
        92.136175603687079,
        95.719694542143202,
        99.330612454787428,
        102.968198614513810,
        106.631760260643450,
        110.320639714757390,
        114.034211781461690,
        117.771881399745060,
        121.533081515438640,
        125.317271149356880,
        129.123933639127240,
        132.952575035616290,
        136.802722637326350,
        140.673923648234250,
        144.565743946344900,
        148.477766951773020,
        152.409592584497350,
        156.360836303078800,
        160.331128216630930,
        164.320112263195170,
        168.327445448427650,
        172.352797139162820,
        176.395848406997370,
        180.456291417543780,
        184.533828861449510,
        188.628173423671600,
        192.739047287844900,
        196.866181672889980,
        201.009316399281570,
        205.168199482641200,
        209.342586752536820,
        213.532241494563270,
        217.736934113954250,
        221.956441819130360,
        226.190548323727570,
        230.439043565776930,
        234.701723442818260,
        238.978389561834350,
        243.268849002982730,
        247.572914096186910,
        251.890402209723190,
        256.221135550009480,
        260.564940971863220,
        264.921649798552780,
        269.291097651019810,
        273.673124285693690,
        278.067573440366120,
        282.474292687630400,
        286.893133295426990,
        291.323950094270290,
        295.766601350760600,
        300.220948647014100,
        304.686856765668720,
        309.164193580146900,
        313.652829949878990,
        318.152639620209300,
        322.663499126726210,
        327.185287703775200,
        331.717887196928470,
        336.261181979198450,
        340.815058870798960,
        345.379407062266860,
        349.954118040770250,
        354.539085519440790,
        359.134205369575340,
        363.739375555563470,
        368.354496072404690,
        372.979468885689020,
        377.614197873918670,
        382.258588773060010,
        386.912549123217560,
        391.575988217329610,
        396.248817051791490,
        400.930948278915760,
        405.622296161144900,
        410.322776526937280,
        415.032306728249580,
        419.750805599544780,
        424.478193418257090,
        429.214391866651570,
        433.959323995014870,
        438.712914186121170,
        443.475088120918940,
        448.245772745384610,
        453.024896238496130,
        457.812387981278110,
        462.608178526874890,
        467.412199571608080,
        472.224383926980520,
        477.044665492585580,
        481.872979229887900,
        486.709261136839360,
        491.553448223298010,
        496.405478487217580,
        501.265290891579240,
        506.132825342034830,
        511.008022665236070,
        515.890824587822520,
        520.781173716044240,
        525.679013515995050,
        530.584288294433580,
        535.496943180169520,
        540.416924105997740,
        545.344177791154950,
        550.278651724285620,
        555.220294146894960,
        560.169054037273100,
        565.124881094874350,
        570.087725725134190,
        575.057539024710200,
        580.034272767130800,
        585.017879388839220,
        590.008311975617860,
        595.005524249382010,
        600.009470555327430,
        605.020105849423770,
        610.037385686238740,
        615.061266207084940,
        620.091704128477430,
        625.128656730891070,
        630.172081847810200,
        635.221937855059760,
        640.278183660408100,
        645.340778693435030,
        650.409682895655240,
        655.484856710889060,
        660.566261075873510,
        665.653857411105950,
        670.747607611912710,
        675.847474039736880,
        680.953419513637530,
        686.065407301994010,
        691.183401114410800,
        696.307365093814040,
        701.437263808737160,
        706.573062245787470,
        711.714725802289990,
        716.862220279103440,
        722.015511873601330,
        727.174567172815840,
        732.339353146739310,
        737.509837141777440,
        742.685986874351220,
        747.867770424643370,
        753.055156230484160,
        758.248113081374300,
        763.446610112640200,
        768.650616799717000,
        773.860102952558460,
        779.075038710167410,
        784.295394535245690,
        789.521141208958970,
        794.752249825813460,
        799.988691788643450,
        805.230438803703120,
        810.477462875863580,
        815.729736303910160,
        820.987231675937890,
        826.249921864842800,
        831.517780023906310,
        836.790779582469900,
        842.068894241700490,
        847.352097970438420,
        852.640365001133090,
        857.933669825857460,
        863.231987192405430,
        868.535292100464630,
        873.843559797865740,
        879.156765776907600,
        884.474885770751830,
        889.797895749890240,
        895.125771918679900,
        900.458490711945270,
        905.796028791646340,
        911.138363043611210,
        916.485470574328820,
        921.837328707804890,
        927.193914982476710,
        932.555207148186240,
        937.921183163208070,
        943.291821191335660,
        948.667099599019820,
        954.046996952560450,
        959.431492015349480,
        964.820563745165940,
        970.214191291518320,
        975.612353993036210,
        981.015031374908400,
        986.422203146368590,
        991.833849198223450,
        997.249949600427840,
        1002.670484599700300,
        1008.095434617181700,
        1013.524780246136200,
        1018.958502249690200,
        1024.396581558613400,
        1029.838999269135500,
        1035.285736640801600,
        1040.736775094367400,
        1046.192096209724900,
        1051.651681723869200,
        1057.115513528895000,
        1062.583573670030100,
        1068.055844343701400,
        1073.532307895632800,
        1079.012946818975000,
        1084.497743752465600,
        1089.986681478622400,
        1095.479742921962700,
        1100.976911147256000,
        1106.478169357800900,
        1111.983500893733000,
        1117.492889230361000,
        1123.006317976526100,
        1128.523770872990800,
        1134.045231790853000,
        1139.570684729984800,
        1145.100113817496100,
        1150.633503306223700,
        1156.170837573242400,
    };
    if (n < 0) { sx_check(SX_UNEXPECTED_VALUE, "Factorial of a negative integer"); return EXIT_FAILURE; }
    else if (n > 254) { x = n + 1; return (x - 0.5) * log(x) - x + 0.5 * log(2 * M_PI) + 1.0 / (12.0 * x); }
    else return lf[n];
}

static double chebyshev_eval(double x, const double *a, const int n)
{
    double b0, b1, b2, twox;
    int i;
    if ((n < 1) || (n > 1000) || (x < -1.1) || (x > 1.1)) return NAN;
    twox = x * 2;
    b2 = b1 = 0;
    b0 = 0;
    for (i = 1; i <= n; i++) { b2 = b1; b1 = b0; b0 = twox * b1 - b2 + a[n - i]; }
    return (b0 - b2) * 0.5;
}

double sx_biomcmc_log1p(double x)
{
    static const double alnrcs[43] = {
        +.10378693562743769800686267719098e+1, -.13364301504908918098766041553133e+0,
        +.19408249135520563357926199374750e-1, -.30107551127535777690376537776592e-2,
        +.48694614797154850090456366509137e-3, -.81054881893175356066809943008622e-4,
        +.13778847799559524782938251496059e-4, -.23802210894358970251369992914935e-5,
        +.41640416213865183476391859901989e-6, -.73595828378075994984266837031998e-7,
        +.13117611876241674949152294345011e-7, -.23546709317742425136696092330175e-8,
        +.42522773276034997775638052962567e-9, -.77190894134840796826108107493300e-10,
        +.14075746481359069909215356472191e-10, -.25769072058024680627537078627584e-11,
        +.47342406666294421849154395005938e-12, -.87249012674742641745301263292675e-13,
        +.16124614902740551465739833119115e-13, -.29875652015665773006710792416815e-14,
        +.55480701209082887983041321697279e-15, -.10324619158271569595141333961932e-15,
        +.19250239203049851177878503244868e-16, -.35955073465265150011189707844266e-17,
        +.67264542537876857892194574226773e-18, -.12602624168735219252082425637546e-18,
        +.23644884408606210044916158955519e-19, -.44419377050807936898878389179733e-20,
        +.83546594464034259016241293994666e-21, -.15731559416479562574899253521066e-21,
        +.29653128740247422686154369706666e-22, -.55949583481815947292156013226666e-23,
        +.10566354268835681048187284138666e-23, -.19972483680670204548314999466666e-24,
        +.37782977818839361421049855999999e-25, -.71531586889081740345038165333333e-26,
        +.13552488463674213646502024533333e-26, -.25694673048487567430079829333333e-27,
        +.48747756066216949076459519999999e-28, -.92542112530849715321132373333333e-29,
        +.17578597841760239233269760000000e-29, -.33410026677731010351377066666666e-30,
        +.63533936180236187354180266666666e-31, };
    if (x == 0.) return 0.;
    if (x == -1.) return INFINITY;
    if (x < -1) return NAN;
    if (fabs(x) <= .375) {
        if (fabs(x) < .5 * DBL_EPSILON) return x;
        if (((0 < x) && (x < 1e-8)) || (-1e-9 < x && x < 0)) return x * (1 - .5 * x);
        return x * (1 - x * chebyshev_eval(x / .375, alnrcs, 22));
    }
    return gsl_log1p(x);
}

double sx_logscale_add(double logx, double logy)
{
    if (logx > logy) return logx + sx_biomcmc_log1p(exp(logy - logx));
    else return logy + sx_biomcmc_log1p(exp(logx - logy));
}

double sx_kahan_logscale_add(double logx, double logy, double *comp_sum)
{
    double corr = 0, nsum = 0, sum = 0, input = 0;
    if (logx > logy) {
        sum = logx;
        input = sx_biomcmc_log1p(exp(logy - logx));
        corr = input - *comp_sum;
        nsum = sum + corr;
        *comp_sum = (nsum - sum) - corr;
        return nsum;
    } else {
        sum = logy;
        input = sx_biomcmc_log1p(exp(logx - logy));
        corr = input - *comp_sum;
        nsum = sum + corr;
        *comp_sum = (nsum - sum) - corr;
        return nsum;
    }
}

double sx_logscale_sub(double logx, double logy)
{
    return logx + sx_biomcmc_log1p(-exp(logy - logx));
}

void sx_skahan_sum(double *sum, double input, double *compensation)
{
    double corr = 0, nsum = 0;
    corr = input - *compensation;
    nsum = *sum + corr;
    *compensation = (nsum - *sum) - corr;
    *sum = nsum;
}

double sx_vkahan_sum(const double *array, int n_elements)
{
    double comp = 0, corr = 0, sum = 0, nsum = 0;
    for (int i = 0; i < n_elements; ++i) {
        corr = array[i] - comp;
        nsum = sum + corr;
        comp = (nsum - sum) - corr;
        sum = nsum;
    }
    return sum;
}

size_t sx_sample_ndoubles(size_t n, const double *array, gsl_rng *seed)
{
    size_t i = 0;
    double sum = 0, value = 0;
    for (i = 0; i < n; ++i) sum += array[i];
    value = gsl_rng_uniform_pos(seed) * sum;
    sum = 0;
    for (i = 0; i < n; ++i) { sum += array[i]; if (sum >= value) break; }
    return i;
}

double sx_kahan_logscale_prob_coal(int i_lin, int o_lin, double bound_time, int Ne)
{
    double logC = 0, logsum_even = 0, logsum_odd = 0, comp_even = 0, comp_odd = 0;
    int i = 0, pi = 0;
    for (i = 0; i < o_lin; ++i) logC += log((o_lin + i) * (i_lin - i) / (double)(i_lin + i));
    logsum_even = (-o_lin * (o_lin - 1) * bound_time / 2 / Ne) + logC;
    switch (i_lin - o_lin) {
    case 0:
        return exp(logsum_even - sx_logscale_fact(o_lin));
    default:
        logC += log((2 * o_lin) * (i_lin - o_lin) / (double)(i_lin + o_lin));
        logsum_odd = (-(o_lin + 1) * o_lin * bound_time / 2 / Ne) + sx_biomcmc_log1p(0.5 / o_lin) + logC;
        for (i = o_lin + 2; i <= i_lin; ++i) {
            pi = i - 1;
            logC += log((o_lin + pi) * (i_lin - pi) / (double)(i_lin + pi) / (i - o_lin));
            if ((i - o_lin) % 2 == 0)
                logsum_even = sx_kahan_logscale_add(logsum_even, (-i * pi * bound_time / 2 / Ne) + 0.693147180559945309417232121458 + log((i - 0.5) / (pi + o_lin)) + logC, &comp_even);
            else
                logsum_odd = sx_kahan_logscale_add(logsum_odd, (-i * pi * bound_time / 2 / Ne) + 0.693147180559945309417232121458 + log((i - 0.5) / (pi + o_lin)) + logC, &comp_odd);
        }
        break;
    }
    if (logsum_odd >= logsum_even) return 0;
    else return exp(sx_logscale_sub(logsum_even, logsum_odd) - sx_logscale_fact(o_lin));
}

/* ---------------- modified Brent (SimPhy brent_root) ---------------- */
typedef double (*root_fn)(double, void *);

static long brent_root(root_fn f, void *ctx, double a, double b, float epsilon, double *result, int verbosity)
{
    int i = 0;
    double prev_b, bprev_b = 0, new_val;
    double fa, fb, fprev_b, fnew_val, initial_a = a, initial_b = b;
    double cond1;
    int prev_bis = 1;
    (void)verbosity;

    fa = f(a, ctx);
    fb = f(b, ctx);
    if (fa * fb >= 0) {
        fprintf(stderr, "Root-finding function error, root (0) not bracketed\n");
        *result = b;
        return SX_UNEXPECTED_VALUE;
    }
    if (fabs(fa) < fabs(fb)) { prev_b = a; a = b; b = prev_b; fprev_b = fa; fa = fb; fb = fprev_b; }
    prev_b = a;
    fprev_b = fa;
    while (((fb != 0) && (fabs(a - b) > epsilon)) || (b == initial_a || b == initial_b)) {
        if ((fa != fprev_b) && (fb != fprev_b))
            new_val = a * fb * fprev_b / (fa - fb) / (fa - fprev_b) + b * fa * fprev_b / (fb - fa) / (fb - fprev_b) + prev_b * fa * fb / (fprev_b - fa) / (fprev_b - fb);
        else
            new_val = b - fb * (b - a) / (fb - fa);
        cond1 = (3 * a + b) / 4;
        if ((!(((new_val > cond1) && (new_val < b)) || ((new_val < cond1) && (new_val > b)))) || (prev_bis && (fabs(new_val - b) >= (fabs(b - prev_b) / 2))) || (!prev_bis && (fabs(new_val - b) >= (fabs(prev_b - bprev_b) / 2)))) {
            new_val = (a + b) / 2;
            prev_bis = 1;
        } else {
            if ((prev_bis && (fabs(b - prev_b) < epsilon)) || (!prev_bis && (fabs(prev_b - bprev_b) < epsilon))) {
                new_val = (a + b) / 2;
                prev_bis = 1;
            } else
                prev_bis = 0;
        }
        fnew_val = f(new_val, ctx);
        bprev_b = prev_b;
        prev_b = b;
        if (fa * fnew_val < 0) { b = new_val; fb = fnew_val; }
        else { a = new_val; fa = fnew_val; }
        if (fabs(fa) < fabs(fb)) { prev_b = a; a = b; b = prev_b; fprev_b = fa; fa = fb; fb = fprev_b; }
        ++i;
        if (i > MAX_IT) { *result = b; return SX_UNEXPECTED_VALUE; }
    }
    *result = b;
    return SX_NO_ERROR;
}

/* ---------------- double-precision attempt ---------------- */
typedef struct {
    int i_lin, o_lin, pop_size;
    double bound_time, lambda_i, lambda_o, C, out_sum, density;
} cdf_ctx;

static double sample_cdf_double(double w_time, void *vctx)
{
    cdf_ctx *cx = vctx;
    double C = cx->C, sum = 0, lambda_k = 0, bound_time = cx->bound_time, lambda_i = cx->lambda_i, lambda_o = cx->lambda_o, out_sum = cx->out_sum, density = cx->density;
    int k = 0, k1 = 0, i_lin = cx->i_lin, o_lin = cx->o_lin, pop_size = cx->pop_size;
    k1 = o_lin - 1;
    if (w_time <= 0) return w_time - density;
    else if (w_time >= bound_time) return 1 - density + (w_time - bound_time);
    else {
        sum = exp(lambda_o * bound_time) * (exp((lambda_i - lambda_o) * w_time) - 1) / (lambda_i - lambda_o) * C;
        for (k = o_lin + 1; k < i_lin; ++k) {
            ++k1;
            lambda_k = -k * k1 / (double)(2 * pop_size);
            C = (o_lin + k1) * (i_lin - 1 - k1) / (double)(i_lin - 1 + k1) / (o_lin - k) * C;
            sum += exp(lambda_k * bound_time) * (exp((lambda_i - lambda_k) * w_time) - 1) / (lambda_i - lambda_k) * (2 * k - 1) / (k1 + o_lin) * C;
        }
        return sum * out_sum - density;
    }
}

double sx_sample_coal_time_mlc(int i_lin, int o_lin, double bound_time, int pop_size, double brent_epsilon, double density, int verbosity)
{
    double lambda_i = 0, lambda_o = 0, C0 = 1, out_sum = 0, ret_value = 0;
    int i = 0;
    lambda_i = -i_lin * (i_lin - 1) / (double)(2 * pop_size);
    lambda_o = -o_lin * (o_lin - 1) / (double)(2 * pop_size);
    for (i = 0; i < o_lin; ++i) C0 *= (o_lin + i) * (i_lin - i - 1) / (double)(i_lin - 1 + i);
    if (o_lin > GSL_SF_FACT_NMAX) return bound_time + 1;
    out_sum = (-lambda_i) / (gsl_sf_fact(o_lin) * sx_kahan_logscale_prob_coal(i_lin, o_lin, bound_time, pop_size));
    cdf_ctx cx = { i_lin, o_lin, pop_size, bound_time, lambda_i, lambda_o, C0, out_sum, density };
    if (brent_root(sample_cdf_double, &cx, 0.0, bound_time, (float)brent_epsilon, &ret_value, verbosity) != SX_NO_ERROR)
        return bound_time + 1;
    else
        return ret_value;
}

/* ---------------- MPFR fallback ---------------- */
static double mpfr_original_prob_coal(int i_lin, int o_lin, double bound_time, int Ne, int precision)
{
    mpfr_t sum, c, c_mult, n_sum, fact;
    int i = 0, j = 0;
    double sum_ret = 0;
    mpfr_inits2(precision, sum, c, c_mult, n_sum, fact, (mpfr_ptr)0);
    mpfr_set_si(sum, 0, MPFR_RNDN);
    for (i = o_lin; i <= i_lin; ++i) {
        mpfr_set_si(c, 1, MPFR_RNDN);
        for (j = 0; j < i; ++j) {
            mpfr_set_si(c_mult, (o_lin + j) * (i_lin - j), MPFR_RNDN);
            mpfr_div_si(c_mult, c_mult, i_lin + j, MPFR_RNDN);
            mpfr_mul(c, c, c_mult, MPFR_RNDN);
        }
        mpfr_set_si(n_sum, -i * (i - 1), MPFR_RNDN);
        mpfr_mul_d(n_sum, n_sum, bound_time, MPFR_RNDN);
        mpfr_div_si(n_sum, n_sum, 2 * Ne, MPFR_RNDN);
        mpfr_exp(n_sum, n_sum, MPFR_RNDN);
        mpfr_mul_si(n_sum, n_sum, 2 * i - 1, MPFR_RNDN);
        mpfr_mul(n_sum, n_sum, c, MPFR_RNDN);
        mpfr_fac_ui(fact, (unsigned long)(i - o_lin), MPFR_RNDN);
        mpfr_div(n_sum, n_sum, fact, MPFR_RNDN);
        mpfr_fac_ui(fact, (unsigned long)o_lin, MPFR_RNDN);
        mpfr_div(n_sum, n_sum, fact, MPFR_RNDN);
        mpfr_div_si(n_sum, n_sum, i + o_lin - 1, MPFR_RNDN);
        if ((i - o_lin) % 2 == 0) mpfr_add(sum, sum, n_sum, MPFR_RNDN);
        else mpfr_sub(sum, sum, n_sum, MPFR_RNDN);
    }
    sum_ret = mpfr_get_ld(sum, MPFR_RNDN);
    mpfr_clears(sum, c, c_mult, n_sum, fact, (mpfr_ptr)0);
    return sum_ret;
}

int sx_mpfr_rescue = 0;
unsigned long sx_mpfr_rescue_count = 0;

typedef struct {
    int i_lin, o_lin, pop_size, precision;
    int out_of_precision;
    double bound_time, density;
    mpfr_ptr lambda_i, lambda_o, initC, out_sum;
    /* temporaries (allocated once per chain) */
    mpfr_t lambda_k, sum, dif_lambdas, nsum, t2, Cx, C;
} mpfr_ctx;

static double sample_cdf_mpfr(double w_time, void *vctx)
{
    mpfr_ctx *cx = vctx;
    double bound_time = cx->bound_time, density = cx->density, cprob = 0;
    int k = 0, k1 = 0, i_lin = cx->i_lin, o_lin = cx->o_lin, pop_size = cx->pop_size;
    k1 = o_lin - 1;
    if (w_time <= 0) return w_time - density;
    else if (w_time >= bound_time) return 1 - density + (w_time - bound_time);
    else {
        mpfr_set(cx->C, cx->initC, MPFR_RNDN);
        mpfr_mul_d(cx->sum, cx->lambda_o, bound_time, MPFR_RNDN);
        mpfr_exp(cx->sum, cx->sum, MPFR_RNDN);
        mpfr_mul(cx->sum, cx->sum, cx->C, MPFR_RNDN);
        mpfr_sub(cx->dif_lambdas, cx->lambda_i, cx->lambda_o, MPFR_RNDN);
        mpfr_div(cx->sum, cx->sum, cx->dif_lambdas, MPFR_RNDN);
        mpfr_mul_d(cx->t2, cx->dif_lambdas, w_time, MPFR_RNDN);
        mpfr_exp(cx->t2, cx->t2, MPFR_RNDN);
        mpfr_sub_si(cx->t2, cx->t2, 1, MPFR_RNDN);
        mpfr_mul(cx->sum, cx->sum, cx->t2, MPFR_RNDN);
        for (k = o_lin + 1; k < i_lin; ++k) {
            ++k1;
            mpfr_set_si(cx->lambda_k, -k * k1, MPFR_RNDN);
            mpfr_div_si(cx->lambda_k, cx->lambda_k, 2 * pop_size, MPFR_RNDN);
            mpfr_set_si(cx->Cx, (o_lin + k1) * (i_lin - 1 - k1), MPFR_RNDN);
            mpfr_div_si(cx->Cx, cx->Cx, i_lin - 1 + k1, MPFR_RNDN);
            mpfr_div_si(cx->Cx, cx->Cx, o_lin - k, MPFR_RNDN);
            mpfr_mul(cx->C, cx->C, cx->Cx, MPFR_RNDN);
            mpfr_mul_d(cx->nsum, cx->lambda_k, bound_time, MPFR_RNDN);
            mpfr_exp(cx->nsum, cx->nsum, MPFR_RNDN);
            mpfr_sub(cx->dif_lambdas, cx->lambda_i, cx->lambda_k, MPFR_RNDN);
            mpfr_mul_d(cx->t2, cx->dif_lambdas, w_time, MPFR_RNDN);
            mpfr_exp(cx->t2, cx->t2, MPFR_RNDN);
            mpfr_sub_si(cx->t2, cx->t2, 1, MPFR_RNDN);
            mpfr_mul(cx->nsum, cx->nsum, cx->t2, MPFR_RNDN);
            mpfr_div(cx->nsum, cx->nsum, cx->dif_lambdas, MPFR_RNDN);
            mpfr_mul_si(cx->nsum, cx->nsum, 2 * k - 1, MPFR_RNDN);
            mpfr_div_si(cx->nsum, cx->nsum, k1 + o_lin, MPFR_RNDN);
            mpfr_mul(cx->nsum, cx->nsum, cx->C, MPFR_RNDN);
            mpfr_add(cx->sum, cx->sum, cx->nsum, MPFR_RNDN);
        }
        mpfr_mul(cx->sum, cx->sum, cx->out_sum, MPFR_RNDN);
        cprob = mpfr_get_d(cx->sum, MPFR_RNDN);
        if (cprob > 1 || cprob < 0) {
            if (!sx_mpfr_rescue) sx_check(SX_UNEXPECTED_VALUE, "MPFR out of precision!!!! That options must be completely out of sense\n");
            cx->out_of_precision = 1;
        }
        return cprob - density;
    }
}

static double mpfr_sample_once(int i_lin, int o_lin, double bound_time, int pop_size, double brent_epsilon, double density, int verbosity, int precision, int *out_of_precision)
{
    mpfr_t lambda_i, lambda_o, C0, Cx, out_sum, fact;
    int i = 0;
    double result = 0;
    mpfr_inits2(precision, lambda_i, lambda_o, C0, Cx, out_sum, fact, (mpfr_ptr)0);
    mpfr_set_si(lambda_i, -i_lin * (i_lin - 1), MPFR_RNDN);
    mpfr_set_si(lambda_o, -o_lin * (o_lin - 1), MPFR_RNDN);
    mpfr_div_si(lambda_i, lambda_i, 2 * pop_size, MPFR_RNDN);
    mpfr_div_si(lambda_o, lambda_o, 2 * pop_size, MPFR_RNDN);
    mpfr_set_si(C0, 1, MPFR_RNDN);
    for (i = 0; i < o_lin; ++i) {
        mpfr_set_si(Cx, (o_lin + i) * (i_lin - i - 1), MPFR_RNDN);
        mpfr_div_si(Cx, Cx, i_lin - 1 + i, MPFR_RNDN);
        mpfr_mul(C0, C0, Cx, MPFR_RNDN);
    }
    mpfr_set(out_sum, lambda_i, MPFR_RNDN);
    mpfr_div_d(out_sum, out_sum, mpfr_original_prob_coal(i_lin, o_lin, bound_time, pop_size, precision), MPFR_RNDN);
    mpfr_fac_ui(fact, (unsigned long)o_lin, MPFR_RNDN);
    mpfr_div(out_sum, out_sum, fact, MPFR_RNDN);
    mpfr_mul_si(out_sum, out_sum, -1, MPFR_RNDN);

    mpfr_ctx cx;
    cx.i_lin = i_lin; cx.o_lin = o_lin; cx.pop_size = pop_size; cx.precision = precision; cx.out_of_precision = 0;
    cx.bound_time = bound_time; cx.density = density;
    cx.lambda_i = lambda_i; cx.lambda_o = lambda_o; cx.initC = C0; cx.out_sum = out_sum;
    mpfr_inits2(precision, cx.lambda_k, cx.sum, cx.dif_lambdas, cx.t2, cx.Cx, cx.nsum, cx.C, (mpfr_ptr)0);
    long rc = brent_root(sample_cdf_mpfr, &cx, 0.0, bound_time, (float)brent_epsilon, &result, verbosity);
    if (!cx.out_of_precision) sx_check(rc, "Error rooting the MPFRCDFCoalTimeFromXtoYLineages\n");
    *out_of_precision = cx.out_of_precision;
    mpfr_clears(cx.lambda_k, cx.sum, cx.dif_lambdas, cx.t2, cx.Cx, cx.nsum, cx.C, (mpfr_ptr)0);
    mpfr_clears(Cx, fact, lambda_i, lambda_o, C0, out_sum, (mpfr_ptr)0);
    return result;
}

double sx_mpfr_sample_coal_time_mlc(int i_lin, int o_lin, double bound_time, int pop_size, double brent_epsilon, double density, int verbosity, int precision)
{
    int oop = 0;
    double r = mpfr_sample_once(i_lin, o_lin, bound_time, pop_size, brent_epsilon, density, verbosity, precision, &oop);
    if (!oop) return r;
    /* rescue: SimPhy 1.0.2 would have aborted here; retry with more precision */
    for (int prec = precision * 2; prec <= 65536; prec *= 2) {
        __sync_fetch_and_add(&sx_mpfr_rescue_count, 1UL);
        r = mpfr_sample_once(i_lin, o_lin, bound_time, pop_size, brent_epsilon, density, verbosity, prec, &oop);
        if (!oop) return r;
    }
    sx_check(SX_UNEXPECTED_VALUE, "MPFR out of precision even at 65536 bits\n");
    return r;
}

void sx_mpfr_chain_run(sx_mpfr_chain_job *job)
{
    double current_ngen = job->n_gen;
    int avail_leaves = job->n_ilin;
    for (int j = 0; j < job->n_coals; ++j) {
        job->out[j] = sx_mpfr_sample_coal_time_mlc(avail_leaves, job->n_olin, current_ngen - job->min_ngen, job->pop_size, job->epsilon, job->u[j], job->verbosity, 512);
        current_ngen -= job->out[j];
        --avail_leaves;
    }
    mpfr_free_cache();
}

double sx_debug_mpfr_cdf(int i_lin, int o_lin, double bound_time, int pop_size, double w, int precision)
{
    mpfr_t lambda_i, lambda_o, C0, Cx, out_sum, fact;
    int i;
    mpfr_inits2(precision, lambda_i, lambda_o, C0, Cx, out_sum, fact, (mpfr_ptr)0);
    mpfr_set_si(lambda_i, -i_lin * (i_lin - 1), MPFR_RNDN);
    mpfr_set_si(lambda_o, -o_lin * (o_lin - 1), MPFR_RNDN);
    mpfr_div_si(lambda_i, lambda_i, 2 * pop_size, MPFR_RNDN);
    mpfr_div_si(lambda_o, lambda_o, 2 * pop_size, MPFR_RNDN);
    mpfr_set_si(C0, 1, MPFR_RNDN);
    for (i = 0; i < o_lin; ++i) {
        mpfr_set_si(Cx, (o_lin + i) * (i_lin - i - 1), MPFR_RNDN);
        mpfr_div_si(Cx, Cx, i_lin - 1 + i, MPFR_RNDN);
        mpfr_mul(C0, C0, Cx, MPFR_RNDN);
    }
    mpfr_set(out_sum, lambda_i, MPFR_RNDN);
    mpfr_div_d(out_sum, out_sum, mpfr_original_prob_coal(i_lin, o_lin, bound_time, pop_size, precision), MPFR_RNDN);
    mpfr_fac_ui(fact, (unsigned long)o_lin, MPFR_RNDN);
    mpfr_div(out_sum, out_sum, fact, MPFR_RNDN);
    mpfr_mul_si(out_sum, out_sum, -1, MPFR_RNDN);
    mpfr_ctx cx;
    cx.i_lin = i_lin; cx.o_lin = o_lin; cx.pop_size = pop_size; cx.precision = precision;
    cx.bound_time = bound_time; cx.density = 0.0;
    cx.lambda_i = lambda_i; cx.lambda_o = lambda_o; cx.initC = C0; cx.out_sum = out_sum;
    mpfr_inits2(precision, cx.lambda_k, cx.sum, cx.dif_lambdas, cx.t2, cx.Cx, cx.nsum, cx.C, (mpfr_ptr)0);
    double v = sample_cdf_mpfr(w, &cx);
    mpfr_clears(cx.lambda_k, cx.sum, cx.dif_lambdas, cx.t2, cx.Cx, cx.nsum, cx.C, (mpfr_ptr)0);
    mpfr_clears(Cx, fact, lambda_i, lambda_o, C0, out_sum, (mpfr_ptr)0);
    return v;
}

double sx_debug_double_cdf(int i_lin, int o_lin, double bound_time, int pop_size, double w)
{
    double lambda_i = -i_lin * (i_lin - 1) / (double)(2 * pop_size);
    double lambda_o = -o_lin * (o_lin - 1) / (double)(2 * pop_size);
    double C0 = 1;
    for (int i = 0; i < o_lin; ++i) C0 *= (o_lin + i) * (i_lin - i - 1) / (double)(i_lin - 1 + i);
    double out_sum = (-lambda_i) / (gsl_sf_fact(o_lin) * sx_kahan_logscale_prob_coal(i_lin, o_lin, bound_time, pop_size));
    cdf_ctx cx = { i_lin, o_lin, pop_size, bound_time, lambda_i, lambda_o, C0, out_sum, 0.0 };
    return sample_cdf_double(w, &cx);
}
