using MultiPrecision;

namespace GaussKronrodTest {
    internal struct Plus32<N> : IConstant where N : struct, IConstant {
        public readonly int Value => checked(default(N).Value + 32);
    }
}
