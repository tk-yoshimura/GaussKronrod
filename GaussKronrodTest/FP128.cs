using MultiPrecision;

namespace GaussKronrodTest {
    internal static class FP128 {
        public static string ToFP128<N>(MultiPrecision<N> x) where N : struct, IConstant {
            Sign sign = x.Sign;
            long exponent = x.Exponent;
            uint[] mantissa = x.Mantissa.Reverse().ToArray();

            string code = $"({(sign == Sign.Plus ? "+1" : "-1")}, {exponent}, 0x{mantissa[0]:X8}{mantissa[1]:X8}uL, 0x{mantissa[2]:X8}{mantissa[3]:X8}uL)";

            return code;
        }
    }
}
