/* 
    ! Editing from GNU 2.7 !

    Original Author:
        Original FORTRAN77 version by Robert Piessens, Maria Branders.
        C++ version by John Burkardt.
    Reference:
        Robert Piessens, Maria Branders,
        A Note on the Optimal Addition of Abscissas to Quadrature Formulas
        of Gauss and Lobatto,
        Mathematics of Computation,
        Volume 28, Number 125, January 1974, pages 135-139.
 */

using MultiPrecision;

namespace GaussKronrod {
    public static class CoefGenaratorMP<N> where N : struct, IConstant {
        static void AbscissaWeights1(int n, int m, MultiPrecision<N> coef2, bool even, MultiPrecision<N>[] b, ref MultiPrecision<N> x, ref MultiPrecision<N> w) {
            MultiPrecision<N> ai, b0 = 0, b1, b2, d0, d1, d2 = 0, delta, dif, f, fd = 0, yy;
            int i, iter, k, ka;

            ka = (x == 0) ? 1 : 0;

            for (iter = 1; iter <= 1024; iter++) {
                b1 = 0;
                b2 = b[m];
                yy = 4 * x * x - 2;
                d1 = 0;

                if (even) {
                    ai = m + m + 1;
                    d2 = ai * b[m];
                    dif = 2;
                }
                else {
                    ai = m + 1;
                    d2 = 0;
                    dif = 1;
                }

                for (k = 1; k <= m; k++) {
                    ai -= dif;
                    i = m - k + 1;
                    b0 = b1;
                    b1 = b2;
                    d0 = d1;
                    d1 = d2;
                    b2 = yy * b1 - b0 + b[i - 1];
                    if (!even) {
                        i++;
                    }
                    d2 = yy * d1 - d0 + ai * b[i - 1];
                }

                if (even) {
                    f = x * (b2 - b1);
                    fd = d2 + d1;
                }
                else {
                    f = 0.5 * (b2 - b0);
                    fd = 4 * x * d2;
                }
                //
                //  Newton correction.
                //
                delta = f / fd;
                x -= delta;

                if (ka == 1) {
                    break;
                }
                if (x.Exponent - delta.Exponent > MultiPrecision<N>.Bits - 8) {
                    ka = 1;
                }
            }

            //
            //  Catch non-convergence.
            //
            if (ka != 1) {
                throw new ArithmeticException("not convergenced");
            }

            //
            //  Computation of the weight.
            //
            d0 = 1;
            d1 = x;
            ai = 0;
            for (k = 2; k <= n; k++) {
                ai += 1;
                d2 = ((ai + ai + 1) * x * d1 - ai * d0) / (ai + 1);
                d0 = d1;
                d1 = d2;
            }

            w = coef2 / (fd * d2);
        }

        static void AbscissaWeights2(int n, int m, MultiPrecision<N> coef2, bool even, MultiPrecision<N>[] b, ref MultiPrecision<N> x, ref MultiPrecision<N> w1, ref MultiPrecision<N> w2) {
            MultiPrecision<N> ai, an, delta, p0 = 0, p1, p2 = 0, pd0, pd1, pd2 = 0, yy;
            int i, iter, k, ka;

            ka = (x == 0) ? 1 : 0;

            //
            //  Iterative process for the computation of a Gaussian abscissa.
            //
            for (iter = 1; iter <= 1024; iter++) {
                p0 = 1;
                p1 = x;
                pd0 = 0;
                pd1 = 1;
                //
                //  When N is 1, we need to initialize P2 and PD2 to avoid problems with DELTA.
                //
                if (n <= 1) {
                    if (MultiPrecision<N>.Epsilon < MultiPrecision<N>.Abs(x)) {
                        p2 = (3 * x * x - 1) / 2;
                        pd2 = 3 * x;
                    }
                    else {
                        p2 = 3 * x;
                        pd2 = 3;
                    }
                }

                ai = 0;
                for (k = 2; k <= n; k++) {
                    ai += 1;
                    p2 = ((ai + ai + 1) * x * p1 - ai * p0) / (ai + 1);
                    pd2 = ((ai + ai + 1) * (p1 + x * pd1) - ai * pd0) / (ai + 1);
                    p0 = p1;
                    p1 = p2;
                    pd0 = pd1;
                    pd1 = pd2;
                }
                //
                //  Newton correction.
                //
                delta = p2 / pd2;
                x -= delta;

                if (ka == 1) {
                    break;
                }
                if (x.Exponent - delta.Exponent > MultiPrecision<N>.Bits - 8) {
                    ka = 1;
                }
            }
            //
            //  Catch non-convergence.
            //
            if (ka != 1) {
                throw new ArithmeticException("not convergenced");
            }
            //
            //  Computation of the weight.
            //
            an = n;

            w2 = 2 / (an * pd2 * p0);

            p1 = 0;
            p2 = b[m];
            yy = 4 * x * x - 2;
            for (k = 1; k <= m; k++) {
                i = m - k + 1;
                p0 = p1;
                p1 = p2;
                p2 = yy * p1 - p0 + b[i - 1];
            }

            if (even) {
                w1 = w2 + coef2 / (pd2 * x * (p2 - p1));
            }
            else {
                w1 = w2 + 2 * coef2 / (pd2 * (p2 - p0));
            }
        }

        public static (MultiPrecision<N>[] x, MultiPrecision<N>[] w1, MultiPrecision<N>[] w2) Coef(int n) {
            MultiPrecision<N>[] x = new MultiPrecision<N>[n + 1], w1 = new MultiPrecision<N>[n + 1], w2 = new MultiPrecision<N>[n + 1];

            MultiPrecision<N> ak, an, bb, c, coef, coef2, d, s, x1, xx, y;

            int i, k, l, ll, m;

            bool even;

            MultiPrecision<N>[] b = new MultiPrecision<N>[((n + 1) / 2) + 1];
            MultiPrecision<N>[] tau = new MultiPrecision<N>[(n + 1) / 2];

            m = (n + 1) / 2;
            even = (2 * m == n);

            d = 2;
            an = 0;
            for (k = 1; k <= n; k++) {
                an += 1;
                d = d * an / (an + 0.5);
            }

            //
            //  Calculation of the Chebyshev coefficients of the orthogonal polynomial.
            //
            tau[0] = (an + 2) / (an + an + 3);
            b[m - 1] = tau[0] - 1;
            ak = an;

            for (l = 1; l < m; l++) {
                ak += 2;
                tau[l] = ((ak - 1) * ak - an * (an + 1)) * (ak + 2) * tau[l - 1]
                    / (ak * ((ak + 3) * (ak + 2) - an * (an + 1)));

                b[m - l - 1] = tau[l];

                for (ll = 1; ll <= l; ll++) {
                    b[m - l - 1] = b[m - l - 1] + tau[ll - 1] * b[m - l + ll - 1];
                }
            }

            b[m] = 1;
            //
            //  Calculation of approximate values for the abscissas.
            //
            bb = MultiPrecision<N>.SinPI(1 / (2 * (an + an + 1)));
            x1 = MultiPrecision<N>.Sqrt(1 - bb * bb);
            s = 2 * bb * x1;
            c = MultiPrecision<N>.Sqrt(1 - s * s);
            coef = 1 - (1 - 1 / an) / (8 * an * an);
            xx = coef * x1;
            //
            //  Coefficient needed for weights.
            //
            //  COEF2 = 2^(2*n+1) * n! * n! / (2n+1)! 
            //        = 2 * 4^n * n! / product( (n+1)*...*(2*n+1))
            //
            coef2 = MultiPrecision<N>.Div(2, 2 * n + 1);
            for (i = 1; i <= n; i++) {
                coef2 = coef2 * 4 * i / (n + i);
            }
            //
            //  Calculation of the K-th abscissa (a Kronrod abscissa) and the
            //  corresponding weight.
            //
            for (k = 1; k <= n; k += 2) {
                AbscissaWeights1(n, m, coef2, even, b, ref xx, ref w1[k - 1]);
                w2[k - 1] = 0;

                x[k - 1] = xx;
                y = x1;
                x1 = y * c - bb * s;
                bb = y * s + bb * c;

                if (k == n) {
                    xx = 0;
                }
                else {
                    xx = coef * x1;
                }

                //
                //  Calculation of the K+1 abscissa (a Gaussian abscissa) and the
                //  corresponding weights.
                //
                AbscissaWeights2(n, m, coef2, even, b, ref xx, ref w1[k], ref w2[k]);

                x[k] = xx;
                y = x1;
                x1 = y * c - bb * s;
                bb = y * s + bb * c;
                xx = coef * x1;
            }
            //
            //  If N is even, we have one more Kronrod abscissa to compute,
            //  namely the origin.
            //
            if (even) {
                xx = 0;
                AbscissaWeights1(n, m, coef2, even, b, ref xx, ref w1[n]);
                w2[n] = 0;
                x[n] = xx;
            }

            return (x, w1, w2);
        }
    }
}
