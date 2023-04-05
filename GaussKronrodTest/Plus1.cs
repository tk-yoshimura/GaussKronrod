using MultiPrecision;

namespace GaussKronrodTest {
    internal struct Plus2<N> : IConstant where N : struct, IConstant {
        public int Value => checked(default(N).Value + 2);
    }
}
